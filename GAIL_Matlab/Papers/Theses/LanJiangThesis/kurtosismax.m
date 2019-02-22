function kmax = kurtosismax(nsig,alphasig,fudge)
kmax =(nsig-3)/(nsig-1) + ((alphasig*nsig)/(1-alphasig))*(1-1/fudge^2)^2;
end