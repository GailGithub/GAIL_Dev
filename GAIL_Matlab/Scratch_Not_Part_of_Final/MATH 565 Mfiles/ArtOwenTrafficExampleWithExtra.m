%Art Owen's Implementation of the Nagel-Schreckenberg Traffic Example
%% Initialization
format compact %eliminate blank lines in output
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20) %make font labels large enough
close all, clear all, %garbage collection

tstart=tic; %start timer
Nvec=20:20:500;
nN=length(Nvec);
totaldist=zeros(nN,1);
for ii=1:nN
    N=Nvec(ii); %number of vehicles
    M=1000; %number of spaces
    vmax=5; %speed limit
    p=1/4; %probability of slowing
    T0=2500; %number of burn in time steps
    T1=1000; %number of time steps to be observed
    T=T0+T1; %total number of time steps

    x=zeros(T+1,N); %initialize vehicle locations
    v=zeros(1,N); %initialize velocity
    x0=randperm(M); %initial placement of vehicles
    x(1,:)=sort(x0(1:N));

    %% Time stepping updates of position
    sumv=0;
    vall=zeros(T+1,N);
    for i=1:T
        d=[diff(x(i,:)) x(i,1)+M-x(i,N)]; %update distances between vehicles
        v=min(v+1,vmax); %speed up if below the speed limit
        v=min(v,d-1); %but do not bump into the vehicle in front
        slowdown=rand(1,N)<p; %perhaps slow down
        v(slowdown)=max(0,v(slowdown)-1);
        x(i+1,:)=x(i,:)+v; %update position of vehicles
        if x(i+1,N) > M %fix last vehicle's position
            x(i+1,N)=x(i+1,N)-M; %last vehicle goes to the start
            x(i+1,:)=[x(i+1,N) x(i+1,1:N-1)]; %shift vector of positions
            v=[v(N) v(1:N-1)];
        end
        sumv=sumv+sum(v);
    end
    display(['Total distance traveled by ' num2str(N) ' cars = ' num2str(sumv)])
    display(['Time for simulation = ' num2str(toc(tstart)) ' seconds'])
    totaldist(ii)=sumv;
end

%% Display results
tstart=tic; %start timer
figure
plot(Nvec,totaldist,'b.','markersize',20)
xlabel('Number of Cars')
ylabel('Total Distance Traveled')
title('Nagel-Schreckenberg Traffic')
display(['Time for the graphical display = ' num2str(toc(tstart)) ' seconds'])


    
    


