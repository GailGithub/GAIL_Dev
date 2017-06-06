function check = funappx_g_check(fappx,out_param)
a = out_param.a;
b = out_param.b;
n = min(ceil(log10(out_param.npoints))+1, 7);
x = a:10^-n:b;
f = out_param.f;
err = abs(f(x)-fappx(x));
[~,ind] = find(err > out_param.abstol);
figure; plot(x, f(x));
figure; plot(x, fappx(x));
figure; semilogy(x, err); xlabel('x'); ylabel('absolute error'); hold on
semilogy(x(ind), err(ind), 'r'); hold off
axis tight
err = max(err);
check = err < out_param.abstol;