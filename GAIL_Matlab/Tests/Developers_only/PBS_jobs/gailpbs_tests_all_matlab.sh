#!/bin/sh
# Driver script to call each matlab version specifically

# delete the TestResults from previous run
rm /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_tests/Test_Results-R*.txt

sh /home/gail/GAIL_tests/PBS_jobs/gitpull.sh

sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2017a"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2017b"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2018a"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2018b"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2019a"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2019b"
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_tests.sh "R2020a"

