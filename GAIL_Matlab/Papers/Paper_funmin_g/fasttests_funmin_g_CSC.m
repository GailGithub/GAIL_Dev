% fasttests_funmin_g_CSC
format
doctest funmin_g_CSC
doctest sinetest
doctest quadratic
doctest flat_bottom
doctest flatbottom2
doctest fejer_jackson_inequality

results=run(ut_funmin_g_CSC)

results=run(ut_workout_funmin_g_CSC)    % time-consuming

run(ut_workout_funmin_g_CSC)   
% Running ut_workout_funmin_g_CSC
%  
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    11  1.20%->55.30%   46.30%      3.70%   44.60%    5.40% 
%   101 32.30%->73.60%   65.90%      3.20%   26.40%    4.50% 
%  1001 66.50%->92.70%   81.90%      4.80%    7.30%    6.00% 
% . 
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    11  1.20%->21.20%   21.20%      0.00%   78.80%    0.00% 
%   101 34.20%->52.40%   52.40%      0.00%   47.60%    0.00% 
%  1001 66.00%->85.30%   80.90%      4.40%   14.70%    0.00% 
% . 
%         Probability    Success   Success   Failure  Failure
%  tau      In Cone    No Warning  Warning No Warning Warning
%    11  2.20%->19.40%   93.60%      0.00%    6.40%    0.00% 
%   101 34.30%->53.60%   95.70%      0.00%    4.30%    0.00% 
%  1001 67.80%->85.30%   96.90%      0.00%    3.10%    0.00% 
%
%            Success    Success    Success    Success
%  TolX               No Warning   Warning    fminbnd
% 1e-02     100.00%    100.00%      0.00%     68.30% 
% 1e-04     100.00%    100.00%      0.00%     68.30% 
% 1e-07     100.00%      0.00%    100.00%     68.30% 