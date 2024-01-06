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

%Lee Frame
EA = 10^6;
EI = 2*10^5;
GA_= 10^5;

D = [EA  0  0;
    0  GA_ 0;
    0   0  EI];

Lee_Frame = twoD_Beam('Lee Frame',10,10,1);
Lee_Frame = set_Constitutive_Matrix_GeoExact(Lee_Frame,D);
Lee_Frame = set_dirichLetBC(Lee_Frame,1);


Lee_Frame_Geometrically_Exact = solveNonLinearEQ(Lee_Frame, 'Geometrically Exact', -133, 'Arch-Length'); % After Bucking
% Lee_Frame_Geometrically_Exact = solveNonLinearEQ(Lee_Frame, 'Geometrically vExact', -132, 'Arch-Length'); % Berfore Bucking
Lee_Frame_Geometrically_Moderate_Rotation = solveNonLinearEQ(Lee_Frame, 'Geometrically Moderate Rotation', -65, 'Arch-Length');
Lee_Frame_Bernoulli = solveNonLinearEQ(Lee_Frame, 'Bernoulli', -150, 'Arch-Length');
Lee_Frame_2nd_order = solveNonLinearEQ(Lee_Frame, '2nd_Order', -17.22, 'Newton-Raphson');


figure(1)
plot_deformedShape(Lee_Frame_Geometrically_Exact, 1, "Lee Frame",1)
plot_lamda_deflection(Lee_Frame_Geometrically_Exact,2)

plot_deformedShape(Lee_Frame_Geometrically_Moderate_Rotation, 1, "Lee Frame",3)
plot_lamda_deflection(Lee_Frame_Geometrically_Moderate_Rotation,4)

plot_deformedShape(Lee_Frame_Bernoulli, 1, "Lee Frame",5)
plot_lamda_deflection(Lee_Frame_Bernoulli,6)


plot_deformedShape(Lee_Frame_2nd_order, 1, "Lee Frame",7)
plot_lamda_deflection(Lee_Frame_2nd_order,8)


%% Plot all the lambda-deflection to geather
figure(2)
store1 = Lee_Frame_Geometrically_Exact.lamda_deflection;
p1 = plot(store1(:,1),store1(:,2),'r-');
set(p1,'DisplayName','Geometrically Exact')

hold on

store2 = Lee_Frame_Geometrically_Moderate_Rotation.lamda_deflection;
p2 = plot(store2(:,1),store2(:,2),'b-');
set(p2, 'DisplayName','Moderate Rotation');
hold on

store3 = Lee_Frame_Bernoulli.lamda_deflection;
p3 = plot(store3(:,1),store3(:,2),'g-');
set(p3,'DisplayName', 'Bernoulli');
hold on


store4 = Lee_Frame_2nd_order.lamda_deflection;
p4 = plot(store4(:,1),store4(:,2),'k-');
set(p4, 'DisplayName','2nd order');

legend
set(gca,'Xdir','reverse')
grid on
axis square
xlabel('Deflection')
ylabel('Load')



opt.FileName = 'plotAxisLimit.jpg'; 
opt.YGrid = 'on'; % 'on' or 'off'
opt.XGrid = 'on'; % 'on' or 'off'
opt.Title = 'Load deflection curve of Lee frame';
opt.XLabel = 'Deflection';
opt.YLabel = 'Load';
% create the plot
setPlotProp(opt);