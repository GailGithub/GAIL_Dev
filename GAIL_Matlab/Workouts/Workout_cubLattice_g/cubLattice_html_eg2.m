% cubLattice_html_eg2.m
format compact
f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
count = 0;
for i=1:100
    [q,out_param] = cubLattice_g(f,hyperbox,'normal',1e-3,1e-3,'transform','C1sin','shift',2^(-25));
    exactsol = 1; check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max');
    if check==0 || isfinite(q) ==0 || out_param.exitflag > 0,
        i, exactsol, q, exitflag = out_param.exitflag,
        abserr = abs(exactsol-q), tol = gail.tolfun(1e-3,1e-3,1,exactsol,'max'), n = out_param.n,
        shift = out_param.shift, lattice = mod(gail.lattice_gen(1,2^24,3)+shift,1);
        max_lattice = max(max(lattice))
        max_C1sin = max(max(lattice-sin(2*pi*lattice)/(2*pi))),
        max_after_normtransform = max(max(gail.stdnorminv(lattice-sin(2*pi*lattice)/(2*pi))))%, min(min(gail.stdnorminv(lattice-sin(2*pi*lattice)/(2*pi))))
          disp('-----');
        count = count + 1;
        %keyboard
        clear lattice
    end;
end;
count

