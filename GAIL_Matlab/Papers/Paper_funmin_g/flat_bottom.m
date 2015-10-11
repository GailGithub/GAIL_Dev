%Widen  output interval of flat bottom function
function flat_bottom;
format compact
format long
f=@(x) exp(-1./(x-.5).^2);
[fmin,out_min] = funmin_g_CSC(f,0,1);
funmin_g_demo(fmin,out_min)
intervals = out_min.intervals % gives [ 0.234036344744525 0.765963655255475 ]
x = out_min.x;

%k = 9
idx1 = find(x<=out_min.intervals(1));
k = length(find(abs(f(x(idx1)) - fmin) < out_min.abstol)) - 1;
xx1 = x(idx1(end-k:end)) % These points xx are not in out_min.intervals
yy1 = f(xx1)             
errors1 = prod(abs(yy1 - fmin) < out_min.abstol) % but abs(f(xx) - fmin) < abstol

%k = 10
idx2 = find(x>=out_min.intervals(2));
k = length(find(abs(f(x(idx2)) - fmin) < out_min.abstol));
xx2 = x(idx2(1:k)) % These points xx are not in out_min.intervals
yy2 = f(xx2)             
errors2 = prod(abs(yy2 - fmin) < out_min.abstol) % but abs(f(xx) - fmin) < abstol

% Shouldn't out_min.intervals be  0.231060606060606   0.768939393939394
[xx1(1), xx2(end)] 