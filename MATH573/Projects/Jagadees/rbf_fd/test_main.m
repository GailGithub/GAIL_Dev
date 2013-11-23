
clear
format compact

debug = false;

u1  = @(x,y) (sin(pi*x).*sin(pi*y));
Lu1 = @(x,y) (-2.0*(pi^2)*sin(pi*x).*sin(pi*y));

u2  = @(x,y) exp(-(x-0.1).^2-0.5*y.^2);
Lu2 = @(x,y) exp(-(x-0.1).^2-0.5*y.^2).*(y.^2 + (-2*x+0.2).^2 - 3);

u3  = @(x,y) (exp(x).*cos(y));
Lu3 = @(x,y) (0);

u4  = @(r,p) (r.^2.*(r-1).*sin(2*p));
Lu4 = @(r,p) (5*r.*sin(2*p));

u5  = @(x,y) (sin(2*x.*y));
Lu5 = @(x,y) (-4*sin(2*x.*y).*((x.^2)+(y.^2)));

u6  = @(x,y) (sin(2*pi*(x-y)));
Lu6 = @(x,y) (-8.0*(pi^2)*sin(2*pi*(x-y)));

u7  = @(x,y) (sin((x.^3).*y) + exp(x) - x./(1+y.^2));
Lu7 = @(x,y) (-9*sin((x.^3).*y).*(x.^4).*(y.^2) + 6*cos((x.^3).*y).*x.*y + exp(x) - sin((x.^3).*y).*(x.^6) - (8*x.*(y.^2))./((1+y.^2).^3) + (2*x)./((1+y.^2).^2)  );

u8  = @(x,y) 0.5*(u1(x,y) + u2(x,y));
Lu8 = @(x,y) 0.5*(Lu1(x,y) + Lu2(x,y));

testU  = u1;
testLu = Lu1;
RBFQR_flag = true;
    
g = 'squareg';
%g = 'circleg';

[p1,p2,p3,p4,p5, pi1,pb1,pi2,pb2,pi3,pb3,pi4,pb4,pi5,pb5, dm1,dm2,dm3,dm4,dm5, np] = gen_tri_mesh(g);

epvec=logspace(log10(0.01),log10(5),50); %debug

err_rbfdir = zeros(length(epvec),5);
stable_flag = zeros(length(epvec),5);
rmsed = zeros(length(epvec),5);

i = 1;
tstart1 = tic;
for epsilon=epvec
    
    if debug==true, fprintf('epsilon %f \n', epsilon); end
    
    [err_rbfdir(i,1), rmsed(i,1), stable_flag(i,1)] = calc_rmse(dm1, p1,pi1,pb1, testU, testU, testLu, epsilon, RBFQR_flag);
    [err_rbfdir(i,2), rmsed(i,2), stable_flag(i,2)] = calc_rmse(dm2, p2,pi2,pb2, testU, testU, testLu, epsilon, RBFQR_flag);
    %[err_rbfdir(i,3), rmsed(i,3), stable_flag(i,3)] = calc_rmse(dm3, p3,pi3,pb3, testU, testU, testLu, epsilon, RBFQR_flag);
    %[err_rbfdir(i,4), rmsed(i,4), stable_flag(i,4)] = calc_rmse(dm4, p4,pi4,pb4, testU, testU, testLu, epsilon, RBFQR_flag);
    
    if debug==true, fprintf('epsilon %1.2f, rmse %1.4f, rmsed %1.4f\n', epsilon, rmse, rmsed), end

    i = i+1;
end
toc(tstart1)

[c,Ieopt1]=min(err_rbfdir(:,1));
[c,Ieopt2]=min(err_rbfdir(:,2));
[c,Ieopt3]=min(err_rbfdir(:,3));
[c,Ieopt4]=min(err_rbfdir(:,4));
[c,Ieopt5]=min(err_rbfdir(:,5));

figure()
loglog(epvec,err_rbfdir)
axis tight
title(sprintf('%s, %s',g,func2str(testU)))
xlabel('\epsilon')
ylabel('RMS error')
ylim([10^-5 100]);
legend(sprintf('E1[%di,%db],\\epsilon%1.2f', np(1,1), np(1,2), epvec(Ieopt1)), ...
       sprintf('E2[%di,%db],\\epsilon%1.2f', np(2,1), np(2,2), epvec(Ieopt2)), ...
       sprintf('E3[%di,%db],\\epsilon%1.2f', np(3,1), np(3,2), epvec(Ieopt3)), ...
       sprintf('E4[%di,%db],\\epsilon%1.2f', np(4,1), np(4,2), epvec(Ieopt4))  );
hold on
x = find(stable_flag(:,1)==true, 1, 'last'); plot(epvec(x), err_rbfdir(x,1), 'r*');
x = find(stable_flag(:,2)==true, 1, 'last'); plot(epvec(x), err_rbfdir(x,2), 'r*');

figure()
loglog(epvec,rmsed)
axis tight
title(sprintf('%s, %s',g,func2str(testU)))
xlabel('\epsilon')
ylabel('Diff RMS error')
legend(sprintf('E1[%di,%db]', np(1,1), np(1,2)), ...
       sprintf('E2[%di,%db]', np(2,1), np(2,2)), ...
       sprintf('E3[%di,%db]', np(3,1), np(3,2)), ...
       sprintf('E4[%di,%db]', np(4,1), np(4,2))  );
hold on
x = find(stable_flag(:,1)==true, 1, 'last'); plot(epvec(x), rmsed(x,1), 'r*');
x = find(stable_flag(:,2)==true, 1, 'last'); plot(epvec(x), rmsed(x,2), 'r*');



