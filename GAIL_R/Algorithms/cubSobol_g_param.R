##Defining the Parameters
cubSobol_g_param = function(r_lag,hyperbox,f,measure, 
                            abstol,reltol,mmin,mmax,fudge,toltype,theta){ 
##To define the variables but we still need to code for the checks****
  #default.hyperbox = [zeros(1,1);ones(1,1)];% default hyperbox
  #default.measure  = 'uniform';
  #default.abstol  = 1e-4;
  #default.reltol  = 1e-2;
  #default.mmin  = 10;
  #default.mmax  = 24;
  #default.fudge = @(m) 5*2.^-m;
  #default.toltype  = 'max';
  #default.theta  = 1;
  f = f;
  hyperbox=hyperbox;
  out_param.measure = measure;
  out_param.abstol = abstol;
  out_param.reltol = reltol;
  out_param.mmin = mmin;
  out_param.mmax = mmax;  
  out_param.fudge = fudge;
  out_param.toltype = toltype;
  out_param.theta = theta;
  out_param.d = size(hyperbox,2);
  out_param=c("measure"=out_param.measure,"abstol"=out_param.abstol,"reltol"=out_param.reltol,"mmin"=out_param.mmin,"mmax"=out_param.mmax,
              "fudge"=out_param.fudge,"toltype"=out_param.toltype,"theta"=out_param.theta, "d"=out_param.d)
  out_list = list("out_param"=out_param, "hyperbox"=hyperbox, "f" = f)
  return(out_list)
}