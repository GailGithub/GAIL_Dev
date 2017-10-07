#!/bin/sh
# driver script to call the python script which sends a combined result for the fasttests

python /home/gail/GAIL_tests/PBS_jobs/sendemail_all_matlab.py "/home/gail/GAIL_tests/PBS_jobs/pbs_reports/gail_daily_tests-$(date +\%Y-\%m-\%d).out"

cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_tests/
# rm Test_Results-R*.txt
