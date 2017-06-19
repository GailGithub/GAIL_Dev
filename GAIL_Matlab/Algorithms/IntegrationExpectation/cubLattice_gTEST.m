% EXAMPLES FROM FILES
g.func = @(x) [10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
g.cv = [8,32/3]; hyperbox= [zeros(1,3);2*ones(1,3)];
q = cubSobol_g(g,hyperbox,'uniform',1e-6,0); exactsol = 128/3;
check = abs(exactsol-q) < 1e-6

g.func = @(x) [10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
g.cv = [8,32/3]; hyperbox= [zeros(1,3);2*ones(1,3)];
w = cubLattice_g(g,hyperbox,'uniform',1e-6,0); exactsol = 128/3;
check = abs(exactsol-w) < 1e-6

f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
y = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-2,'transform','C1'); exactsol = (sqrt(pi)/2*(erf(2)+erf(1)))^2;
check = abs(exactsol-y) < max(1e-3,1e-2*abs(exactsol))

f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
q = cubSobol_g(f,hyperbox,'uniform',1e-3,1e-2); exactsol = (sqrt(pi)/2*(erf(2)+erf(1)))^2;
check = abs(exactsol-q) < max(1e-3,1e-2*abs(exactsol))

% -----------------------------------------------------------------
numIter=100; 

%% Sobol' Sampling regular

AMeanCallSobol = optPrice(AMeanCallIID); %make a copy of the IID optPrice object
AMeanCallSobol.priceParam.cubMethod = 'Sobol'; %change to Sobol sampling
pointsSobol(1,numIter) = 0;
for n=1:numIter
    gail.TakeNote(n,10)
    [AMeanCallSobolPrice,AoutSobol] = genOptPrice(AMeanCallSobol);
    pointsSobol(n) = AMeanCallSobolPrice;
end

%% Sobol' Sampling with Control Variates
% We can use control variates with Sobol' and lattice sampling, but it is a
% bit different than for IID sampling.  Here is an example.

AMeanCallSobolCV = optPrice(AMeanCallSobol); %make a copy of the object
AMeanCallSobolCV.payoffParam = struct( ...
    'optType',{{'amean','gmean'}},...  %Add two payoffs
    'putCallType',{{'call','call'}});  %both calls
AMeanCallSobolCV.priceParam.cubMethod = 'SobolCV'; %change method to use control variates
pointsSobolCV(1,numIter) = 0;
for n=1:numIter
    gail.TakeNote(n,10)
    [AMeanCallSobolCVPrice, AoutSobolCV] = genOptPrice(AMeanCallSobolCV);
    pointsSobolCV(n) = AMeanCallSobolCVPrice;
end

%% Lattice Sampling with Control Variates
% We can use control variates with lattice sampling, but it is a
% bit different than for IID sampling.  Here is an example.

AMeanCallLatticeCV = optPrice(AMeanCallLattice); %make a copy of the object
AMeanCallLatticeCV.payoffParam = struct( ...
    'optType',{{'amean','gmean'}},...  %Add two payoffs
    'putCallType',{{'call','call'}});  %both calls
AMeanCallLatticeCV.priceParam.cubMethod = 'LatticeCV'; %change method to use control variates
[AMeanCallLatticeCVPrice,AoutLatticeCV] = genOptPrice(AMeanCallLatticeCV);

pointsLatticeCV(1,numIter) = 0;
for n= 1:numIter
    gail.TakeNote(n,10)
    [AMeanCallLatticeCVPrice,AoutLatticeCV] = genOptPrice(AMeanCallLatticeCV);    
    pointsLatticeCV(n) = AMeanCallLatticeCVPrice;
end

%% Graphs 
avSobol=mean(pointsSobolCV)
rangeSobol=range(pointsSobolCV)
stdSobol=std(pointsSobolCV)

avLat=mean(pointsLatticeCV)
rangeLat=range(pointsLatticeCV)
stdLat=std(pointsLatticeCV)

ax1 = subplot(2,1,1);
x = linspace(0,numIter);
y1 = pointsSobolCV;
hold on; 
plot(ax1,x,y1, 'k')
plot(ax1, x, avSobol);

ax2 = subplot(2,1,2);
y2 = pointsLatticeCV;
plot(ax2,x,y2, 'b')

%% Stem leaf plot
pointsSobolCV*1000;
pointsLatticeCV*1000;
stemleafplot(pointsSobolCV);
stemleafplot(pointsLatticeCV);