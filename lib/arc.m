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
function [ddl1,ddl2] = arc(da,dab,dat,dl,iq,psi,dll)
       
    c1 = dat'*dat + psi^2*(iq'*iq);
    c2 = 2*(da+dab)'*dat + 2*dl*(psi^2)*(iq'*iq);
    c3 = (da+dab)'*(da+dab) + (dl^2)*(psi^2)*(iq'*iq)-dll^2;
    
    delta = c2^2-4*c1*c3;
    
    if delta>0
        %dls = roots([c1, c2, c3]);
        %ddl1 = dls(1);
        %ddl2 = dls(2);
        ddl1 = (-c2-sqrt(delta))/(2*c1);
        ddl2 = (-c2+sqrt(delta))/(2*c1);
        
        if ddl1>ddl2
            ddl1 = (-c2+sqrt(delta))/(2*c1);
            ddl2 = (-c2-sqrt(delta))/(2*c1);
        end
    else
        ddl1 = -c2 / 2 * c1;
        ddl2 = -c2 / 2 * c1;
        disp('Possible issue in Arc Length equation');
    end
end