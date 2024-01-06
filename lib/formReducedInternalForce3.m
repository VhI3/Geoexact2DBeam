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
function  [reducedInternalForce]=formReducedInternalForce3(GDof,prescribedDof,numberElements,elementNodes,elementAngelS,elementLengthS,D,poly,reducedU)
activeDof = setdiff([1:GDof]',[prescribedDof]);
internalForce = zeros(GDof,1);
U = zeros(GDof,1);
U(activeDof,1) = reducedU;
for e = 1:numberElements
    indice = elementNodes(e,:);
    cos_alpha = cos(elementAngelS(e,:));   sin_alpha = sin(elementAngelS(e,:));    
    ttt = [cos_alpha sin_alpha 0;
          -sin_alpha cos_alpha 0;
            0          0       1];
    T_I = [ttt      zeros(3);
           zeros(3)   ttt];   
    elementDof = [];
    for node = 1:size(indice,2)
        elementDof = [elementDof 3*(indice(node)-1)+1 3*(indice(node)-1)+2 3*(indice(node)-1)+3];
    end
    U_e_l = T_I*U(elementDof);
    Le = elementLengthS(e);
    w_I = U_e_l(2:3:size(U_e_l,1));
    psi_I = U_e_l(3:3:size(U_e_l,1));  
    k_l = zeros(6);      
    wp_I = [w_I(1) psi_I(1) w_I(2) psi_I(2)];
    K11 = zeros(2);
    K22L = zeros(4);
    K22nL = zeros(4);
    K12 = zeros(2,4);
    K21 = zeros(4,2);      
    [ xi , w ] = Gauss_int(poly+1);
    for p = 1:poly+1
        cubic = Cubic(xi(p));
        cubic(2,:) = (2/Le)*cubic(2,:);
        cubic(3,:) = (2/Le)^2*cubic(3,:);        
        cubic(:,2) = Le/2*cubic(:,2);
        cubic(:,4) = Le/2*cubic(:,4);
        N_d = 2/Le*Lagrange1D_d(xi(p),poly);
        w_d_e = 0;
        for JJ = 1:2*size(elementNodes,2)
            w_d_e = w_d_e + cubic(2,JJ)*wp_I(JJ);
        end 
        % Loop for K11, 2X2
        for i = 1:2
            for j = 1:2
                K11(i,j) = K11(i,j) + D(1,1)*N_d(i)*N_d(j)*Le/2*w(p);
            end
        end        
        % Loop for K22(Linear Part), 4X4
        for I = 1:4
            for J = 1:4
                K22L(I,J) = K22L(I,J) + D(3,3)*cubic(3,I)*cubic(3,J)*Le/2*w(p);
            end
        end    
    end          
    [ xi_r , w_r ] = Gauss_int(poly);
    for p = 1:poly
        cubic = Cubic(xi_r(p));
        cubic(2,:) = 2/Le*cubic(2,:);
        cubic(3,:) = (2/Le)^2*cubic(3,:);
        cubic(:,2) = Le/2*cubic(:,2);
        cubic(:,4) = Le/2*cubic(:,4);
        N_d_r = 2/Le*Lagrange1D_d(xi_r(p),poly);
        w_d_e = 0;
        for JJ = 1:2*size(elementNodes,2)
            w_d_e = w_d_e + cubic(2,JJ)*wp_I(JJ);
        end 
         % Loop for K12, 2X4
        for i = 1:2
            for J = 1:4
                K12(i,J) = K12(i,J) + 1/2*D(1,1)*w_d_e*N_d_r(i)*cubic(2,J)*Le/2*w_r(p); 
            end
        end        
         % Loop for K21, 4X2
        for I = 1:4
            for j = 1:2
                K21(I,j) = K21(I,j) + D(1,1)*w_d_e*N_d_r(j)*cubic(2,I)*Le/2*w_r(p);
            end
        end               
        % Loop for K22 and TangentMatrix (Nonlinear part), 4X4
        for I = 1:4
            for J = 1:4
                K22nL(I,J) = K22nL(I,J) + 1/2*D(1,1)*w_d_e^2*cubic(2,I)*cubic(2,J)*Le/2*w_r(p);
            end
        end
    end
    K22 = K22L + K22nL;   
    K = [K11 K12;
         K21 K22];
    Kc(1,:) = K(1,:);
    Kc(2,:) = K(3,:);
    Kc(3,:) = K(4,:);
    Kc(4,:) = K(2,:);
    Kc(5,:) = K(5,:);
    Kc(6,:) = K(6,:);    
    k_l(:,1) = Kc(:,1);
    k_l(:,2) = Kc(:,3);
    k_l(:,3) = Kc(:,4);
    k_l(:,4) = Kc(:,2);
    k_l(:,5) = Kc(:,5);
    k_l(:,6) = Kc(:,6);    
    R_l = k_l*U_e_l;    
    R_G = T_I'*R_l;
    internalForce(elementDof) = internalForce(elementDof) + R_G;
end
reducedInternalForce = internalForce(activeDof);
end