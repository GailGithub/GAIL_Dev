#!/bin/sh

#echo "Tests gave errors. See the output attached." | mutt -a "/home/gail/GAIL_tests/test_results.txt" -s "subject of message" -- gail_dev@googlegroups.com

# attach the file, instead of putting inte body, as the file contains non printable chars 
fileattach="/home/gail/GAIL_tests/PBS_jobs/pbs_reports/gail_daily_tests-$(date +\%Y-\%m-\%d).out"
echo Todays test output file is $fileattach
# Sending PBS reports
uuencode $fileattach $fileattach
mail -s "Daily test results were WRONG: REPORT" -a $fileattach schoi32@iit.edu,jrathin1@iit.edu < /dev/null


# And if you want mail to read the content from a file:
nl=$'\n\n';
begin=" We got some errors in our daily tests:";
uuencode Test_Results.txt Test_Results.txt
# changes=$(< Test_Results.txt);
echo "${begin}${nl}" | mail -s "Daily test results were WRONG: ERRORS" -a "Test_Results.txt" schoi32@iit.edu,jrathin1@iit.edu

#echo "${begin}${nl}${changes}" | mail -s "Daily test results were WRONG: ERRORS" schoi32@iit.edu,jrathin1@iit.edu
#(echo "${begin}${nl}${changes}";uuencode Test_Results.txt Test_Results.txt) | mail -s "GAIL daily test results with ERRORS" gail_dev@googlegroups.com -- -r "ljimene1@hawk.iit.edu"

# Reference
# http://www.binarytides.com/linux-mailx-command/
# https://unix.stackexchange.com/questions/42145/shell-script-hangs-on-mail-command
