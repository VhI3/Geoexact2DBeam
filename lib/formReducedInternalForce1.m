% /*
% MIT License
% 
% Copyright (c) 2024 VHI3
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to use,
% copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
% Software, and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
% THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% */
%%
% Function to form the reduced internal force vector for a beam element.
function [reducedInternalForce] = ...
    formReducedInternalForce1(GDof, prescribedDof, numberElements, ...
    elementNodes, elementAngelS, elementLengthS, D, poly, reducedU)

    % Initialize variables and vectors.
    activeDof = setdiff([1:GDof]', [prescribedDof]);  % Active degrees of freedom (excluding prescribed ones).
    internalForce = zeros(GDof,1);                   % Initialize the internal force vector.
    [xi_r, w_r] = Gauss_int(poly);                    % Gaussian integration points and weights.
    U = zeros(GDof,1);                               % Initialize the displacement vector.
    U(activeDof,1) = reducedU;                        % Assign reduced displacements to active DOFs.

    % Loop over each element to calculate its contribution to the internal force.
    for e = 1:numberElements
        indice = elementNodes(e,:);                   % Nodes associated with the current element.
        cos_alpha = cos(elementAngelS(e,:));          % Cosine of element angle.
        sin_alpha = sin(elementAngelS(e,:));          % Sine of element angle.
        
        % Transformation matrix from local to global coordinates.
        T_I = [cos_alpha sin_alpha 0;
              -sin_alpha cos_alpha 0;
                0          0       1];

        % Initialize element degree of freedom array.
        elementDof = [];
        for node = 1:size(indice,2)
            elementDofInMatrix(node,:) =  [3*(indice(node)-1)+1 3*(indice(node)-1)+2 3*(indice(node)-1)+3];
            elementDof = [elementDof 3*(indice(node)-1)+1 3*(indice(node)-1)+2 3*(indice(node)-1)+3];
            U_l(elementDofInMatrix(node,:),1) = T_I*U(elementDofInMatrix(node,:)',1);
        end
        U_e_l = U_l(elementDof');
         
        % Extract local displacements (u, w) and rotations (psi) for each node in the element.
        u_I = U_e_l(1:3:size(U_e_l,1));
        w_I = U_e_l(2:3:size(U_e_l,1));
        psi_I = U_e_l(3:3:size(U_e_l,1));
        
        % Loop over each node to compute internal forces.
        for I = 1:size(elementNodes,2)
            R_l = zeros(3,1);  % Initialize local force vector.
            for p = 1:poly
                 % Calculate shape function and its derivative.
                 N = Lagrange1D(xi_r(p), poly);
                 N_d = 2/elementLengthS(e) * Lagrange1D_d(xi_r(p), poly);

                 % Compute local displacements and rotations.
                 psi_e = 0; u_d_e = 0; w_d_e = 0;
                 for i = 1:size(elementNodes,2)
                    psi_e = psi_e + N(i)*psi_I(i);
                    u_d_e = u_d_e + N_d(i)*u_I(i);
                    w_d_e = w_d_e + N_d(i)*w_I(i);
                 end
                 c = cos(psi_e);
                 s = sin(psi_e);
                 
                 % Compute strain and stress components.
                 alpha_e = -(1+u_d_e)*s + w_d_e*c;  % Alpha component of strain.
                 betta_e = -(1+u_d_e)*c - w_d_e*s; % Beta component of strain.
                 T_psi = [c s 0;
                         -s c 0;
                          0 0 1]; 
                 NN = [1-c; s; 0];
                 B0_U_I = zeros(3,1);
                 for m = 1:size(elementNodes,2)
                    B0_U_I = B0_U_I + [N_d(m) 0 0;
                                       0 N_d(m) 0;
                                       0 0 N_d(m)] * U_l(elementDofInMatrix(m,:));
                 end
                 epsilon_e = (T_psi*B0_U_I-NN);    % Local strain.
                 S_e = D*epsilon_e;                % Local stress using constitutive matrix D.

                 % Compute the B matrix.
                 B = [N_d(I)*c N_d(I)*s N(I)*alpha_e;
                     -N_d(I)*s N_d(I)*c N(I)*betta_e;
                      0      0       N_d(I)];
                  
                 % Accumulate local force vector.
                 R_l = R_l + (B'*S_e*elementLengthS(e)/2*w_r(p));
            end
            % Transform local forces to global coordinate system.
            R_G = T_I'*R_l;
            internalForce(elementDofInMatrix(I,:)) = internalForce(elementDofInMatrix(I,:)) + R_G;
        end
    end
    % Extract reduced internal force for the active degrees of freedom.
    reducedInternalForce = internalForce(activeDof);
end
