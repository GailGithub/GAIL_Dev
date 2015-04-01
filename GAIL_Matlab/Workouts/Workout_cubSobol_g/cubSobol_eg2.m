format compact
for i=1:100
    q = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3);
    exactsol = 1; check = abs(exactsol-q) < gail.tolfun(1e-3,1e-3,1,exactsol,'max');
    if check==0, 
        i, exactsol, q, abs(exactsol-q), gail.tolfun(1e-3,1e-3,1,exactsol,'max')
        keyboard
    end;
end;