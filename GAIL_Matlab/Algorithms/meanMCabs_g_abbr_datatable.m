meanMCabs_g_abbr()
meanMCabs_g_abbr()
meanMCabs_g_abbr()

ntot_matrix = zeros(10,5)
time_matrix = zeros(10,5)

for n = 1:50
   % if(n<=10)
      %  [tmu,out_param] = meanMCabs_g_abbr(@(n) rand(n,1).^2, 0.0001)
      %  ntot_matrix(n,1) = out_param.ntot
      %  time_matrix(n,1) = out_param.time
    if(n<=20)
        [tmu,out_param] = meanMCabs_g_abbr(@(n) rand(n,1).^2, 0.0005)
        ntot_matrix(n,2) = out_param.ntot
        time_matrix(n,2) = out_param.time

    elseif(n<=30)
        [tmu,out_param] = meanMCabs_g_abbr(@(n) rand(n,1).^2, 0.001)
        ntot_matrix(n,3) = out_param.ntot
        time_matrix(n,3) = out_param.time
    elseif(n<=40)
        [tmu,out_param] = meanMCabs_g_abbr(@(n) rand(n,1).^2, 0.005)
        ntot_matrix(n,4) = out_param.ntot
        time_matrix(n,4) = out_param.time
    elseif(n<=50)
        [tmu,out_param] = meanMCabs_g_abbr(@(n) rand(n,1).^2, 0.01)
        ntot_matrix(n,5) = out_param.ntot
        time_matrix(n,5) = out_param.time
    end
end
format long g
    avg_ntot_matrix = sum(ntot_matrix)/10
    avg_time_matrix = sum(time_matrix)/10
    overall_matrix = [avg_ntot_matrix;avg_time_matrix]