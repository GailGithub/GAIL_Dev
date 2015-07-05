source("meanMCabs_g_abbr.R")
meanMCabs_g_abbr()
meanMCabs_g_abbr()
meanMCabs_g_abbr()

data_0.0001 = c()
data_0.0005 = c()
data_0.001 = c()
data_0.005 = c()
data_0.01 = c()

for (i in 1:50) {
  if (i<=10) { a = 0.0001 
               data_0.0001 = append(data_0.0001,meanMCabs_g_abbr(abstol = a))
  }
    else if(i<=20) { a = 0.0005 
                     data_0.0005 = append(data_0.0005,meanMCabs_g_abbr(abstol = a))
    }
      else if(i<=30) { a = 0.001 
                       data_0.001 = append(data_0.001,meanMCabs_g_abbr(abstol = a))
      }
        else if(i<=40) { a = 0.005 
                         data_0.005 = append(data_0.005,meanMCabs_g_abbr(abstol = a))
        }
          else if(i<=50) { a = 0.01
                           data_0.01 = append(data_0.01,meanMCabs_g_abbr(abstol = a))
          }
}
avg_ntot1 = (data_0.0001[2] + data_0.0001[6] + data_0.0001[10] + data_0.0001[14] + data_0.0001[18] + data_0.0001[22] + data_0.0001[26] + data_0.0001[30] + data_0.0001[34] + data_0.0001[38])/10
avg_time1 = (data_0.0001[4] + data_0.0001[8] + data_0.0001[12] + data_0.0001[16] + data_0.0001[20] + data_0.0001[24] + data_0.0001[28] + data_0.0001[32] + data_0.0001[36] + data_0.0001[40])/10

avg_ntot2 = (data_0.0005[2] + data_0.0005[6] + data_0.0005[10] + data_0.0005[14] + data_0.0005[18] + data_0.0005[22] + data_0.0005[26] + data_0.0005[30] + data_0.0005[34] + data_0.0005[38])/10
avg_time2 = (data_0.0005[4] + data_0.0005[8] + data_0.0005[12] + data_0.0005[16] + data_0.0005[20] + data_0.0005[24] + data_0.0005[28] + data_0.0005[32] + data_0.0005[36] + data_0.0005[40])/10

avg_ntot3 = (data_0.001[2] + data_0.001[6] + data_0.001[10] + data_0.001[14] + data_0.001[18] + data_0.001[22] + data_0.001[26] + data_0.001[30] + data_0.001[34] + data_0.001[38])/10
avg_time3 = (data_0.001[4] + data_0.001[8] + data_0.001[12] + data_0.001[16] + data_0.001[20] + data_0.001[24] + data_0.001[28] + data_0.001[32] + data_0.001[36] + data_0.001[40])/10

avg_ntot4 = (data_0.005[2] + data_0.005[6] + data_0.005[10] + data_0.005[14] + data_0.005[18] + data_0.005[22] + data_0.005[26] + data_0.005[30] + data_0.005[34] + data_0.005[38])/10
avg_time4 = (data_0.005[4] + data_0.005[8] + data_0.005[12] + data_0.005[16] + data_0.005[20] + data_0.005[24] + data_0.005[28] + data_0.005[32] + data_0.005[36] + data_0.005[40])/10

avg_ntot5 = (data_0.01[2] + data_0.01[6] + data_0.01[10] + data_0.01[14] + data_0.01[18] + data_0.01[22] + data_0.01[26] + data_0.01[30] + data_0.01[34] + data_0.01[38])/10
avg_time5 = (data_0.01[4] + data_0.01[8] + data_0.01[12] + data_0.01[16] + data_0.01[20] + data_0.01[24] + data_0.01[28] + data_0.01[32] + data_0.01[36] + data_0.01[40])/10

ntot = c(avg_ntot1,avg_ntot2,avg_ntot3,avg_ntot4,avg_ntot5)
time = c(avg_time1,avg_time2,avg_time3,avg_time4,avg_time5)

overall_data = matrix(c(ntot,time),ncol=5,byrow=TRUE)
colnames(overall_data) = c("0.0001","0.0005","0.001","0.005","0.01")
rownames(overall_data) = c("ntot","time")
overall_data

