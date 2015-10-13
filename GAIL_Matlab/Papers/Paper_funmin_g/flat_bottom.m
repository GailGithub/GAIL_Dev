%Widen  output interval of flat bottom function
function flat_bottom;
format compact
format long
f=@(x) exp(-1./(x-.5).^2);
[fmin,out_min] = funmin_g_CSC(f,0,1);
funmin_g_demo(fmin,out_min)
intervals = out_min.intervals % gives [ 0.234036344744525 0.765963655255475 ]
x = out_min.x;

y = f(x);
[~, index] = find(abs(y - fmin) < out_min.abstol);
leftint = find([1 diff(index)~=1]);
rightint = find([diff(index)~=1 1]);
[x(index(leftint)), x(index(rightint))] 
% Shouldn't out_min.intervals be [0.231060606060606   0.768939393939394]
 