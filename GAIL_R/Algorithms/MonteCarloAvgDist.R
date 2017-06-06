#Given two sets of locations computes the Euclidean distance matrix among pairings.
source('C:/Users/jmate/Desktop/RUnit/meanMC_CLT.R')       #meanMC_CLT.R script
n=1e4;                                              #number of samples.
d=2;                                                #dimensions
print("++++++  OUTPUT ++++++");
     start = Sys.time();                     #TIC.
     storeMeans = vector(mode="numeric", length = 0); #we created a numeric vector of inicial length equals to 0.
     for(i in 1:6){                          #compute the mean in a loop of range equals to i.
                #(X             -        Y)      
                dist=sqrt(rowSums(((matrix(runif(d*n),n,d)-matrix(runif(d*n),n,d))^2)));
                meandist = mean(dist);              #compute the sample average distance of the n rows.
                storeMeans = c(storeMeans, meandist); #store means of each loop in a vector.
                print(meandist);  
      }
      finish = Sys.time() - start;            #TOC.
      print(finish);                          #print elapsed time.
      #print(storeMeans);
      plot(storeMeans,main="meanDist",xlab="quantity", ylab="range") #construct the graph of z.
            
#output using meanMC_CLT function from meanMC_CLT.R
print("++++++ OUTPUT meanMC_CLT ++++++")
out_meanMC_CLT = meanMC_CLT(Yrand = function(n) {sqrt(rowSums(((matrix(runif(d*n),n,d)-matrix(runif(d*n),n,d))^2)))},0.02)
print(out_meanMC_CLT)

#Authors: Marcela Ribeiro, Ramon Oliveira and Joao Mateus Cunha.