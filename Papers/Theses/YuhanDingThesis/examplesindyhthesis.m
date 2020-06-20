function examplesindyhthesis
% To replicate the examples presented in Yuhan Ding's dissertation
% Example 1(the same example in the cones not balls paper)
[succnowarn,succwarn]=conepaper_test_funappx_g(10000,10^7,1e-8);
% From Example 2 to Example 7
% Comparsion with chebfun
chebfuncomparison
% Output the number of points for funappxglobal_g
out_param
% Example to illustrate the computational cost
TableForLocalAlgorithm
% GUI example of a flat function
funappx_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-2,5,9,'funappxPenalty_g');
% Approximate symbol of batman with tolerance 1e-4, nlo as 10, nhi as 15
[negerr, poserr]=newplotbatmanfunappx(1e-4);
max(max(negerr),max(poserr))
% Plot Cyclone Example with default nlo 10 and nhi 100;
a = 0; b = 16*pi; tol = 1e-6; plotcyclone(a,b,tol);
% Approximate a seashell in Florida
a=-0.2; b=0.5; c=0.1; n=2; az=-150; el=10; res=128;
funappxseashell(a, b, c, n, az, el, res);
% Comparison between funappx_g and funappxglobal_g
nrep = 100; abstol = 1e-7; nlo = 100; nhi = 1000;
workout_funappx_g(nrep,abstol,nlo,nhi,'funappxPenalty_g');
