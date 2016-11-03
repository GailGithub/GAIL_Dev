#!/bin/bash
#
# check_GAIL_tests_report.sh
#
#  Created on: Mar 1, 2015
#      Author: Sou-Cheng Choi, sctchoi@uchicago.edu
# 
#
# Run "./check_GAIL_tests_report.sh gail_workouts-2015-02-23-22-00-11.txt" 
#     (without quotes).
# Or run "./check_GAIL_tests_report.sh gail_workouts-2015-02-23-22-00-11.txt 
#       > Test.log 2>&1"  (without quotes) to save 
#      screen output of this shell script to Test.log.
#      Summary of Unit-test run-time errors are found in "Test.log.
# 

cd /home/gail/GAIL_tests/workout_reports

echo "~~~~~~~~~~~~~~~~~~~~~~Summary of errors~~~~~~~~~~~~~~~~~~~~~~"
grep -n Error "$1"  
grep -n "not ok" "$1" 
echo "~~~~~~~~~~~~~~~~~~~~~~~~last few lines~~~~~~~~~~~~~~~~~~~~~~~~"
tail "$1" 
echo "~~~~~~~~~~~~~~~~~~~~total number of lines~~~~~~~~~~~~~~~~~~~~~"
wc -l "$1"
