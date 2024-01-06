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
addpath("lib")
clear variables
% clc
close all

% William Toggle
A = 0.408282;
I = 0.0132651;
E = 19971400;
A_ = 5/6*A;
EA = E*A;
EI = E*I;
nu = 0.27;
G = E/2/(1+nu);
GA_ = G*A_;

D = [EA  0  0;
     0  GA_ 0;
     0   0  EI];

william_toggle = twoD_Beam('William Toggle', 8,8,1);
william_toggle = set_Constitutive_Matrix_GeoExact(william_toggle, D);
william_toggle = set_dirichLetBC(william_toggle,2);

william_toggle_Geometrically_Exact = solveNonLinearEQ(william_toggle,'Geometrically Exact',-1.6,'Arch-Length');
william_toggle_Geometrically_Moderate_Rotation = solveNonLinearEQ(william_toggle, 'Geometrically Moderate Rotation', -1.6,'Arch-Length');
william_toggle_Bernoulli = solveNonLinearEQ(william_toggle, 'Bernoulli', -0.55,'Arch-Length');
william_toggle_2nd_Order = solveNonLinearEQ(william_toggle, '2nd_Order', -0.5,'Newton-Raphson');


plot_deformedShape(william_toggle_Geometrically_Exact, 1, "William Toggle",1)
plot_lamda_deflection(william_toggle_Geometrically_Exact,2)

plot_deformedShape(william_toggle_Geometrically_Moderate_Rotation, 1, "William Toggle",3)
plot_lamda_deflection(william_toggle_Geometrically_Moderate_Rotation,4)

plot_deformedShape(william_toggle_Bernoulli, 1, "William Toggle",5)
plot_lamda_deflection(william_toggle_Bernoulli,6)

plot_deformedShape(william_toggle_2nd_Order, 1, "William Toggle",7)
plot_lamda_deflection(william_toggle_2nd_Order,8)

%% Plot all the lambda-deflection to geather
figure(2)
store1 = william_toggle_Geometrically_Exact.lamda_deflection;
p1 = plot(store1(:,1),store1(:,2),'r-');
set(p1,'DisplayName','Geometrically Exact')

hold on

store2 = william_toggle_Geometrically_Moderate_Rotation.lamda_deflection;
p2 = plot(store2(:,1),store2(:,2),'b-');
set(p2, 'DisplayName','Moderate Rotation');
hold on

store3 = william_toggle_Bernoulli.lamda_deflection;
p3 = plot(store3(:,1),store3(:,2),'g-');
set(p3,'DisplayName', 'Bernoulli');
hold on


store4 = william_toggle_2nd_Order.lamda_deflection;
p4 = plot(store4(:,1),store4(:,2),'k-');
set(p4, 'DisplayName','2nd order');

legend
set(gca,'Xdir','reverse')
grid on
axis square
xlabel('Deflection')
ylabel('Load')