%Plot discrepancy function
clear all, close all

nplot=100;
[xplt,yplt]=meshgrid([0:nplot]/nplot);
xxplt=[xplt(:) yplt(:)];
[nnplt,d]=size(xxplt);

xlatt=mod(([0:7]'*[1 3]+0.5)/8,1); %good lattice
[n,d]=size(xlatt);

Funif=ones(nnplt,1);
temp=ones(nnplt,n);
for j=1:d;
    Funif=Funif.*xxplt(:,j);
    temp=temp.*(repmat(xxplt(:,j),1,n)>=repmat(xlatt(:,j)',nnplt,1));
end
Femp=mean(temp,2);
discfungood=Funif-Femp;
discfunplt=reshape(discfungood,nplot+1,nplot+1);

[mindisc,whmin]=min(discfungood);
[maxdisc,whmax]=max(discfungood);
absmean=mean(abs(discfungood));
rms=sqrt(mean(discfungood.*discfungood));
maxabs=max(abs([maxdisc mindisc]));
disp('Discrepancy function for good lattice takes on values between')
disp(['     ' num2str(mindisc) ' at (x,y) = (' num2str(xxplt(whmin,:)) ') and'])
disp(['     ' num2str(maxdisc) ' at (x,y) = (' num2str(xxplt(whmax,:)) ')'])
disp(['     It has mean absolute value ' num2str(absmean) ])
disp(['         root mean square value ' num2str(rms) ])
disp(['         maximum absolute value ' num2str(maxabs) ])
disp(['                L_2 discrepancy ' num2str(l2disc(xlatt)) ])

figure; h=plot(xlatt(:,1),xlatt(:,2),'rs'); set(h,'linewidth',6);
axis('equal'), axis([0 1 0 1])
xlabel('\it{x_1}')
ylabel('\it{x_2}')

figure; surf(xplt,yplt,discfunplt); 
shading interp; colormap(autumn);

xlatt=mod(([0:7]'*[1 7]+0.5)/8,1); %bad lattice
[n,d]=size(xlatt);

Funif=ones(nnplt,1);
temp=ones(nnplt,n);
for j=1:d;
    Funif=Funif.*xxplt(:,j);
    temp=temp.*(repmat(xxplt(:,j),1,n)>=repmat(xlatt(:,j)',nnplt,1));
end
Femp=mean(temp,2);
discfunbad=Funif-Femp;
discfunplt=reshape(discfunbad,nplot+1,nplot+1);

figure; h=plot(xlatt(:,1),xlatt(:,2),'rs'); set(h,'linewidth',6);
axis('equal'), axis([0 1 0 1])
xlabel('\it{x_1}')
ylabel('\it{x_2}')

figure; surf(xplt,yplt,discfunplt); 
shading interp; colormap(spring);

[mindisc,whmin]=min(discfunbad);
[maxdisc,whmax]=max(discfunbad);
absmean=mean(abs(discfunbad));
rms=sqrt(mean(discfunbad.*discfunbad));
maxabs=max(abs([maxdisc mindisc]));
disp('Discrepancy function for a bad lattice takes on values between')
disp(['     ' num2str(mindisc) ' at (x,y) = (' num2str(xxplt(whmin,:)) ') and'])
disp(['     ' num2str(maxdisc) ' at (x,y) = (' num2str(xxplt(whmax,:)) ')'])
disp(['     It has mean absolute value ' num2str(absmean) ])
disp(['         root mean square value ' num2str(rms) ])
disp(['         maximum absolute value ' num2str(maxabs) ])
disp(['                L_2 discrepancy ' num2str(l2disc(xlatt)) ])

xlatt=rand(8,2); %random numbers
[n,d]=size(xlatt);

Funif=ones(nnplt,1);
temp=ones(nnplt,n);
for j=1:d;
    Funif=Funif.*xxplt(:,j);
    temp=temp.*(repmat(xxplt(:,j),1,n)>=repmat(xlatt(:,j)',nnplt,1));
end
Femp=mean(temp,2);
discfunrnd=Funif-Femp;
discfunplt=reshape(discfunrnd,nplot+1,nplot+1);

figure; h=plot(xlatt(:,1),xlatt(:,2),'rs'); set(h,'linewidth',6);
axis('equal'), axis([0 1 0 1])
xlabel('\it{x_1}')
ylabel('\it{x_2}')

figure; surf(xplt,yplt,discfunplt); 
shading interp; colormap(spring);

[mindisc,whmin]=min(discfunrnd);
[maxdisc,whmax]=max(discfunrnd);
absmean=mean(abs(discfunrnd));
rms=sqrt(mean(discfunrnd.*discfunrnd));
maxabs=max(abs([maxdisc mindisc]));
disp('Discrepancy function for random points takes on values between')
disp(['     ' num2str(mindisc) ' at (x,y) = (' num2str(xxplt(whmin,:)) ') and'])
disp(['     ' num2str(maxdisc) ' at (x,y) = (' num2str(xxplt(whmax,:)) ')'])
disp(['     It has mean absolute value ' num2str(absmean) ])
disp(['         root mean square value ' num2str(rms) ])
disp(['         maximum absolute value ' num2str(maxabs) ])
disp(['                L_2 discrepancy ' num2str(l2disc(xlatt)) ])

figure; h=plot(sort(discfungood),([1:nnplt]-0.5)/nnplt,'-b', sort(discfunbad),([1:nnplt]-0.5)/nnplt,'--k', sort(discfunrnd),([1:nnplt]-0.5)/nnplt,'-.r'); 
set(h,'linewidth',2)
xlabel('Discrepancy Function');
ylabel('Probability');
refresh
