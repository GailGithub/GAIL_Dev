%Art Owen's Implementation of the Nagel-Schreckenberg Traffic Example
%% Initialization
format compact %eliminate blank lines in output
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20) %make font labels large enough
close all, clear all, %garbage collection

tstart=tic; %start timer
N=100; %number of vehicles
M=1000; %number of spaces
vmax=5; %speed limit
p=2/3; %probability of slowing
T0=2500; %number of burn in time steps
T1=1000; %number of time steps to be observed
T=T0+T1; %total number of time steps
flowmax=N*T1*vmax; %maximum flow of cars during the time observed

x=zeros(T+1,N); %initialize vehicle locations
v=zeros(1,N); %initialize velocity to zero
%v=ceil(vmax*rand(1,N)); %initialize velocity randomly
x0=randperm(M); %initial placement of vehicles
x(1,:)=sort(x0(1:N)); %sorted in order

%% Time stepping updates of position
for i=1:T
    d=mod([diff(x(i,:)) x(i,1)-x(i,N)],M); %update distances between vehicles
    v=min(v+1,vmax); %speed up if below the speed limit
    v=min(v,d-1); %but do not bump into the vehicle in front
    slowdown=rand(1,N)<p; %which cars slow down
    v(slowdown)=max(0,v(slowdown)-1); %slow these down
    x(i+1,:)=x(i,:)+v; %update position of vehicles
end
avgvelocity=sum(x(T+1,:)-x(T0+1,:))/(N*T1); %Average velocity of all cars

%% Display results
display(['Time for simulation = ' num2str(toc(tstart)) ' seconds'])
display(['Average velocity of ' num2str(N) ' cars = ' num2str(avgvelocity)])
display(['     which is ' num2str(100*avgvelocity/vmax) ...
    '% of the ideal of ' num2str(vmax)])
display(['Flux = ' num2str(avgvelocity*N) ' cars per unit time'])

tstart=tic; %start timer
figure
plot(mod(x(T0+1:T+1,:),M),repmat((0:T1)',1,N),'b.','markersize',3)
xlabel('Position')
ylabel('Time')
title('Nagel-Schreckenberg Traffic')
print -depsc NSTraffic.eps
display(['Time for the graphical display = ' num2str(toc(tstart)) ' seconds'])


    
    


