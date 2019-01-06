#!/bin/sh
# driver script to call the python script which sends a combined result for the fasttests
source /home/gail/.bashrc
# /home/gail/usr/bin/mega-whoami
# /home/gail/usr/bin/mega-cmd 2>&1 > /home/gail/megacmd.log &
python /home/gail/GAIL_tests/PBS_jobs/sendemail_all_matlab.py "/home/gail/GAIL_tests/PBS_jobs/pbs_reports/gail_daily_tests-$(date +\%Y-\%m-\%d).out"

cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_tests/
# rm Test_Results-R*.txt
