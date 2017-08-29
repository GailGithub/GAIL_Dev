% cubSobol_eg2.m
format compact
f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
count = 0;
for i=1:100
    [q,out_param] = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3,'mmin',21);
    exactsol = 1; check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max');
    if check==0 || isfinite(q) ==0 %|| out_param.exitflag > 0, 
        i, exactsol, q, exitflag = out_param.exitflag, 
        abserr = abs(exactsol-q), tol = gail.tolfun(1e-3,1e-3,1,exactsol,'max')
        disp('-----');
        count = count + 1;
        keyboard
    end;
end;
count