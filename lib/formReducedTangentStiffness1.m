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
function  [reducedTstiffness]=...
    formReducedTangentStiffness1(GDof,prescribedDof,numberElements,...
    elementNodes,elementAngelS,elementLengthS,D,poly,reducedU)

activeDof = setdiff([1:GDof]',[prescribedDof]);
Tstiffness = zeros(GDof); 
[ xi , w ] = Gauss_int(poly);
U = zeros(GDof,1);
U(activeDof,1) = reducedU;
for e = 1:numberElements 
    % elementDof: element degrees of freedom (Dof)
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
        for K = 1:size(elementNodes,2)
            k_l = zeros(3);

           for p = 1:poly
              N = Lagrange1D(xi(p),poly);
              N_d = 2/elementLengthS(e)*Lagrange1D_d(xi(p),poly);

              psi_e = 0;
              u_d_e = 0;
              w_d_e = 0;
              for i = 1:size(elementNodes,2)
                psi_e = psi_e + N(i)*psi_I(i);
                u_d_e = u_d_e + N_d(i)*u_I(i);
                w_d_e = w_d_e + N_d(i)*w_I(i);
              end
             c = cos(psi_e);
             s = sin(psi_e);
             % for alpha_e
             alpha_e = -(1+u_d_e)*s + w_d_e*c;
             % for betta_e
             betta_e = -(1+u_d_e)*c - w_d_e*s;
              % for alpha3
             alpha3 = -(1+u_d_e)*c - w_d_e*s;
             % for alpha4
             alpha4 = (1+u_d_e)*s - w_d_e*c;
             T_psi = [c s 0;
                     -s c 0;
                      0 0 1]; 
             NN = [1-c ;
                    s  ;
                    0 ];
             B0_U_I = zeros(3,1)   ;
             for m = 1: size(elementNodes,2)
                B0_U_I = B0_U_I + [N_d(m)  0   0;
                                    0  N_d(m)  0;
                                    0    0   N_d(m)]*U_l(elementDofInMatrix(m,:));
             end
             S_e = D*(T_psi* B0_U_I - NN);
             N_e = S_e(1);
             Q_e = S_e(2);
             
             B_I = [N_d(I)*c N_d(I)*s N(I)*alpha_e;
                   -N_d(I)*s N_d(I)*c N(I)*betta_e;
                      0          0    N_d(I)     ];
                  
             B_K = [N_d(K)*c N_d(K)*s N(K)*alpha_e;
                   -N_d(K)*s N_d(K)*c N(K)*betta_e;
                      0          0    N_d(K)     ];
                  
             G_N = [     0              0          -N_d(I)*N(K)*s;
                         0              0           N_d(I)*N(K)*c;
                   -N(I)*N_d(K)*s   N(I)*N_d(K)*c   N(I)*N(K)*alpha3];
               
             G_Q = [     0              0          -N_d(I)*N(K)*c;
                         0              0          -N_d(I)*N(K)*s;
                   -N(I)*N_d(K)*c  -N(I)*N_d(K)*s   N(I)*N(K)*alpha4];
             

             k_l = k_l + (B_I'*D*B_K + N_e*G_N + Q_e*G_Q)*elementLengthS(e)/2*w(p);  
           end
           k_G = T_I'*k_l*T_I;
           Row = elementDofInMatrix(I,:);
           Col = elementDofInMatrix(K,:);
           Tstiffness(Row,Col) = Tstiffness(Row,Col) + k_G;
        end
    end
    
end
reducedTstiffness = Tstiffness(activeDof,activeDof);
end
