#echo "Tests gave errors. See the output attached." | mutt -a "/home/gail/GAIL_tests/test_results.txt" -s "subject of message" -- gail_dev@googlegroups.com

# And if you want mail to read the content from a file:
nl=$'\n\n';
begin=" We got some errors in our daily tests:";
changes=$(< Test_Results.txt);
echo "${begin}${nl}${changes}" | mail -s "GAIL daily test results with ERRORS" gail_dev@googlegroups.com -- -r "ljimene1@hawk.iit.edu"

#(echo "${begin}${nl}${changes}";uuencode Test_Results.txt Test_Results.txt) | mail -s "GAIL daily test results with ERRORS" gail_dev@googlegroups.com -- -r "ljimene1@hawk.iit.edu"
