#echo "Tests gave errors. See the output attached." | mutt -a "/home/lantoni/GAIL_tests/test_results.txt" -s "subject of message" -- gail_dev@googlegroups.com

# And if you want mail to read the content from a file:
nl=$'\n\n';
begin=" We got some errors in our daily tests:";
changes=$(< Test_Results.txt);
mail -s "GAIL daily test results with ERRORS" lluisantoni@gmail.com -- -r "ljimene1@hawk.iit.edu" <<< "${begin}${nl}${changes}"
# gail_dev@googlegroups.com, schoi32@iit.edu -- -r "ljimene1@hawk.iit.edu" <<< "${begin}${nl}${changes}"
