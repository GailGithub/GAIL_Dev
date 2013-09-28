% Sandwich simulation

%% Initialize
close all, clear all
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)

%% Parameters
tic;
whole=4; %wholesale price of sandwich
retail=5; %retail price of sandwich
order=20; %quantity of sandwiches ordered daily
demandlo=5; %lo end of demand
demandhi=35; %hi end of demand
ndays=3*365; %number of days for simulation
ndvec=(1:ndays)'; %vector of 1 to number of days
nreps=2500; %number of reptitions
nrvec=(1:nreps)'; %number of reptitions of the simulation

%% Perform simulation
demand=randi([demandlo,demandhi],nreps,ndays); %uniform random numbers for demand
sold=min(demand,order); %amount of sandwiches sold that day
dayprofit=sold*retail-whole*order; %profit for the day
avgdayprofit=cumsum(dayprofit,2)./repmat(ndvec',nreps,1); %average profit for the first n days
toc

%% Output results
%Numerical output
meanday=mean(dayprofit,1)';
meanavg=mean(avgdayprofit,1)';
stderrday=std(dayprofit,0,1)'/sqrt(nreps); %standard error for profit
stderravg=std(avgdayprofit,0,1)'/sqrt(nreps); %standard error for profit
disp(['For ' int2str(nreps) ' replications of'])
disp(['    ' int2str(ndays) ' days of business'])
disp(['For a supply of ' int2str(order) ' sandwiches ordered daily'])
disp('and a random demand that is uniform over a range of')
disp(['   {' int2str(demandlo) ',...,' ...
    int2str(demandhi), '} sandwiches'])
disp(['The average daily profit over this whole time = $' ...
    num2str(meanavg(ndays),' %6.2f') ...
    ' +/- ' num2str(2.58*stderravg(ndays),' %6.2f')])
disp(['   compared to the maximum possible profit of $' ...
    num2str(order*(retail-whole),' %6.2f')])
disp(' ');

%Plot daily and cumulative average profit
subplot(2,1,1)
semilogx(ndvec,meanday,'b-',...
    ndvec,repmat(meanday,1,2)+stderrday*(2.58*[-1 1]),'r--',...
    'linewidth',2)
set(gca,'xtick',10.^(0:ceil(log10(ndays))));
xlabel('Number of Days'); ylabel('Daily Profit')
subplot(2,1,2)
semilogx(ndvec,meanavg,'g',...
    ndvec,repmat(meanavg,1,2)+stderravg*(2.58*[-1 1]),'r--',...
    'linewidth',2);
set(gca,'xtick',10.^(0:ceil(log10(ndays))));
xlabel('Number of Days'); ylabel('Avg Daily Profit')


    