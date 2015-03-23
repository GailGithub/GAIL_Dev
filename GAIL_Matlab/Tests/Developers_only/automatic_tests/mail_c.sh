# echo "Daily test results were OK" | mail -s "GAIL test results" gail_dev@googlegroups.com -- -r "ljimene1@hawk.iit.edu"
mail -s "Daily test results were OK (report)" schoi32@iit.edu, lluisantoni@gmail.com < /home/gail/GAIL_tests/PBS_jobs/pbs_reports/gail_daily_tests-$(date +\%Y-\%m-\%d).out

# (echo "Daily test results were OK";uuencode Test_Results.txt Test_Results.txt) | mail -s "GAIL test results" gail_dev@googlegroups.com -- -r "ljimene1@hawk.iit.edu"
# And if you want mail to read the content from a file:
# mail -s "GAIL daily test results" gail_dev@googlegroups.com </home/gail/GAIL_tests/mail_c.txt
