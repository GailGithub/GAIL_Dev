#!/bin/sh

# attach the file, instead of putting inte body, as the file contains non printable chars 
file_attach="/home/gail/GAIL_tests/PBS_jobs/pbs_reports/gail_daily_tests-$(date +\%Y-\%m-\%d).out"
uuencode $file_attach $file_attach 
mail -s "Daily test results were OK (report)" -a $file_attach schoi32@iit.edu,jrathin1@iit.edu
