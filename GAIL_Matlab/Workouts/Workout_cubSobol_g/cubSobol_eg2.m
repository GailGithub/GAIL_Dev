format compact
f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
for i=1:100
    [q,out_param] = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3);
    exactsol = 1; check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max');
    if check==0, 
        i, exactsol, q, out_param.exitflag, abs(exactsol-q), gail.tolfun(1e-3,1e-3,1,exactsol,'max')
        %keyboard
    end;
end;