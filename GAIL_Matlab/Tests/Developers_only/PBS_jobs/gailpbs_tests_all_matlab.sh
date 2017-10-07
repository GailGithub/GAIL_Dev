#!/bin/sh
# Driver script to call each matlab version specifically

# delete the TestResults from previous run
rm /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_tests/Test_Results-R*.txt

sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2014a"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2014b"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2015a"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2015b"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2016a"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2016b"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2017a"
