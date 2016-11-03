% this is the Nagel-Schreckenberg traffic model 
% input --- replication number/sample size R
% output --- the average velocity of all the cars during the given period of time
function b = Ytrafficmodel(R) %input replication number and output average velocity
vmax = 7;% max velocity
N = 100;% number of cars
M = 1000;% number of spaces
T = 2500;% total number of time steps
T0 = 1500;% the number of time steps that NOT count
p = 1/3;% random stop probability
v = zeros(R,N); % initial velocity for all cars
%v = round(vmax*rand(R,N)); % initial velocity for all cars
x = zeros(R,N);% positions for all cars
for j = 1:R
    RandPosition = randperm(M);
    x(j,:) = sort(RandPosition(1:N))-1; 
    % place N cars into M positions for each replication
end
for t = 1:T
    d = [diff(x,1,2) x(:,1)+M-x(:,N)]; %update distances between vehicles
    v = min(v + 1,vmax); %speed up if below the speed limit
    v = min(v,d-1); %not bump into the vehicle in front
    z = rand(R,N) < p;
    v = max((v-z),0);%random slow down with probability p
    x = x+v;% update the position
    if t == T0
        x0 = x; 
        %If taking T0 as the starting time , then take the x0 as
        %the initial position
    end
end
xT = x;
b = mean(xT-x0,2)/(T-T0);% calculate the average velocity
