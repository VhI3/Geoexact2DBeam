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
function [ Return ] = Cubic(xi)
    Cubic1 = 1/4*(2-3*xi+xi^3); 
    Cubic2 = 1/4*(1-xi-xi^2+xi^3);
    Cubic3 = 1/4*(2+3*xi-xi^3);
    Cubic4 = 1/4*(-1-xi+xi^2+xi^3);
    H = [Cubic1  Cubic2 Cubic3 Cubic4];
    
    Cubic_d1 = 1/4*(-3+3*xi^2); 
    Cubic_d2 = 1/4*(-1-2*xi+3*xi^2);
    Cubic_d3 = 1/4*( 3-3*xi^2);
    Cubic_d4 = 1/4*(-1+2*xi+3*xi^2);
    H_d  = [Cubic_d1 Cubic_d2 Cubic_d3 Cubic_d4];
    
    
    Cubic_dd1 = 1/4*(6*xi);
    Cubic_dd2 = 1/4*(-2+6*xi);
    Cubic_dd3 = 1/4*(-6*xi); 
    Cubic_dd4 = 1/4*( 2+6*xi);
    H_dd = [Cubic_dd1 Cubic_dd2 Cubic_dd3 Cubic_dd4];
        
    
    Return = [H;H_d;H_dd];
end