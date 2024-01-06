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
function [ Nd ] = lagrange_func_d( xi, xik, p )

if(size(xik,2)~=(p+1))          %check the compatible
    return;
end

Nd = zeros(p+1,size(xi,2));      %initial lize the shape functions matrix

for j=1:size(xi,2)%loop over all the xi values
    for i=1:p+1                 %loop over the shape functions
        Nd(i,j)=0;
        for m=1:p+1
            if(m ~= i)
                temp = 1;
                for k=1:p+1                       %loop k
                    if(k ~= i && k ~= m)          %and skip k=i and k=m
                        temp=temp*(xik(k)-xi(j))/(xik(k)-xik(i));
                    end
                end
                Nd(i,j) = Nd(i,j) - temp/(xik(m)-xik(i));
            end
        end
    end
end

end

