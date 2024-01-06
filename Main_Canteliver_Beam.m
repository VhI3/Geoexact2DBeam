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
clc
close all

% Canteliver Beam
EA = 10^6;
EI = 2*10^5;
GA_ = 10^5; 
 D = [EA  0  0;
     0  GA_ 0;
     0   0  EI];


canteliver_beam = twoD_Beam('Canteliver Beam',10,0,1);
canteliver_beam = set_Constitutive_Matrix_GeoExact(canteliver_beam, D);
canteliver_beam = set_dirichLetBC(canteliver_beam,3);


canteliver_beam_Geometrically_Exact = solveNonLinearEQ(canteliver_beam,'Geometrically Exact',-240,'Arch-Length');
canteliver_beam_Geometrically_Moderate_Rotation = solveNonLinearEQ(canteliver_beam,'Geometrically Moderate Rotation',-240,'Arch-Length');
canteliver_beam_Bernoulli = solveNonLinearEQ(canteliver_beam,'Bernoulli',-220,'Arch-Length');
canteliver_beam_2nd_Order = solveNonLinearEQ(canteliver_beam,'2nd_Order',-150,'Newton-Raphson');

plot_deformedShape(canteliver_beam_Geometrically_Exact, 1, "Canteliver Beam",1)
plot_lamda_deflection(canteliver_beam_Geometrically_Exact,2)

plot_deformedShape(canteliver_beam_Geometrically_Moderate_Rotation, 1, "Canteliver Beam",3)
plot_lamda_deflection(canteliver_beam_Geometrically_Moderate_Rotation,4)


plot_deformedShape(canteliver_beam_Bernoulli, 1, "Canteliver Beam",5)
plot_lamda_deflection(canteliver_beam_Bernoulli,6)

plot_deformedShape(canteliver_beam_2nd_Order, 1, "Canteliver Beam",7)
plot_lamda_deflection(canteliver_beam_2nd_Order,8)
%% Plot all the lambda-deflection to geather
figure(2)
store1 = canteliver_beam_Geometrically_Exact.lamda_deflection;
p1 = plot(store1(:,1),store1(:,2),'r-');
set(p1,'DisplayName','Geometrically Exact')

hold on

store2 = canteliver_beam_Geometrically_Moderate_Rotation.lamda_deflection;
p2 = plot(store2(:,1),store2(:,2),'b-');
set(p2, 'DisplayName','Moderate Rotation');
hold on

store3 = canteliver_beam_Bernoulli.lamda_deflection;
p3 = plot(store3(:,1),store3(:,2),'g-');
set(p3,'DisplayName', 'Bernoulli');
hold on


store4 = canteliver_beam_2nd_Order.lamda_deflection;
p4 = plot(store4(:,1),store4(:,2),'k-');
set(p4, 'DisplayName','2nd order');

legend
set(gca,'Xdir','reverse')
grid on
axis square
xlabel('Deflection')
ylabel('Load')


