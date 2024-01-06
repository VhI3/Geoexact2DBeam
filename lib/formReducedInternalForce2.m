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
function  [reducedInternalForce]=...
    formReducedInternalForce2(GDof,prescribedDof,numberElements,...
    elementNodes,elementAngelS,elementLengthS,D,poly,reducedU)
activeDof = setdiff([1:GDof]',[prescribedDof]);
internalForce = zeros(GDof,1);
[ xi_r , w_r ] = Gauss_int(poly);
U = zeros(GDof,1);
U(activeDof,1) = reducedU;
for e = 1:numberElements
    indice = elementNodes(e,:);
    cos_alpha = cos(elementAngelS(e,:));
    sin_alpha = sin(elementAngelS(e,:));    
    T_I = [cos_alpha sin_alpha 0;
          -sin_alpha cos_alpha 0;
            0          0       1];
    elementDof = [];
    for node = 1:size(indice,2)
        elementDofInMatrix(node,:) =  [3*(indice(node)-1)+1 3*(indice(node)-1)+2 3*(indice(node)-1)+3];
        elementDof = [elementDof 3*(indice(node)-1)+1 3*(indice(node)-1)+2 3*(indice(node)-1)+3];
        U_l(elementDofInMatrix(node,:),1) = T_I*U(elementDofInMatrix(node,:)',1);
    end
    U_e_l = U_l(elementDof');
    u_I = U_e_l(1:3:size(U_e_l,1));
    w_I = U_e_l(2:3:size(U_e_l,1));
    psi_I = U_e_l(3:3:size(U_e_l,1));
    for I = 1:size(elementNodes,2)
        R_l = zeros(3,1);
        for p = 1:poly
           N = Lagrange1D(xi_r(p),poly);
           N_d = 2/elementLengthS(e)*Lagrange1D_d(xi_r(p),poly);
           psi_e = 0;
           u_d_e = 0;
           w_d_e = 0;
           for i = 1:size(elementNodes,2)
              psi_e = psi_e + N(i)*psi_I(i);
              u_d_e = u_d_e + N_d(i)*u_I(i);
              w_d_e = w_d_e + N_d(i)*w_I(i);
           end         
           B0_U_I = zeros(3,1);
             for m = 1: size(elementNodes,2)
                B0_U_I = B0_U_I + [N_d(m)  N_d(m)*psi_e  -1/2*psi_e*N(m);
                                    0      N_d(m)        -(1+u_d_e)*N(m);
                                    0       0                N_d(m)     ]*U_l(elementDofInMatrix(m,:));
             end           
           S_e = D*B0_U_I;           
           B = [N_d(I)        N_d(I)*psi_e   (w_d_e-psi_e)*N(I);
               -N_d(I)*psi_e  N_d(I)         -(1+u_d_e)*N(I)   ;
                 0            0                N_d(I)         ];           
           R_l = R_l + (B'*S_e*elementLengthS(e)/2*w_r(p));
        end
        R_G = T_I'*R_l;
        internalForce(elementDofInMatrix(I,:)) = internalForce(elementDofInMatrix(I,:)) + R_G;
    end
end
reducedInternalForce = internalForce(activeDof);
end