function  [pappx1,out_param1,pappx2,out_param2]=parabolatestfunction(tol,nlo,nhi)
%tt=0.1; a=0; b=1; tol=1e-5; nlo=10; nhi=1000;
p1=@(x) -(x-0.5).^2 +25;
p1doubleprime=@(x) -2;
p2=@(x,tt) -5*(x-0.5).^2 +25;
p2doubleprime=@(x,tt) -10;
[pappx1,out_param1]=funappx_g(p1,0,1,tol,nlo,nhi);
[pappx2,out_param2]=funappx_g(p2,0,1,tol,nlo,nhi);
% t = 0:0.0001:1;
% figure;
% subplot(2,1,1);
% plot(t,p1(t));
% title('test function p1(x)')
% subplot(2,1,2);
% plot(t,abs(p1doubleprime(t)));
% title('|p1"(x)|')
% ylim([0,3]);
% gail.save_eps('WorkoutFunappxOutput', 'testp1');
% figure(2);
% subplot(2,1,1);
% plot(t,p2(t));
% title('test function p2(x)')
% subplot(2,1,2);
% plot(t,abs(p2doubleprime(t)));
% title('|p2"(x)|')
% ylim([0,11]);
% gail.save_eps('WorkoutFunappxOutput', 'testp2');
