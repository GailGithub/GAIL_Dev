% cubMC_eg2.m
format compact
f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
count = 0;
for i=1:100
    [q,out_param] = cubMC_g(f,hyperbox,'normal',1e-3,1e-3);
    exactsol = 1; check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max');
    if check==0 || isfinite(q) ==0 || out_param.exit > 1,
        i, exactsol, q, exitflag = out_param.exit,
        abserr = abs(exactsol-q), tol = gail.tolfun(1e-3,1e-3,1,exactsol,'max')
        disp('-----');
        count = count + 1;
        keyboard
    else
        i
    end;
end;
count
