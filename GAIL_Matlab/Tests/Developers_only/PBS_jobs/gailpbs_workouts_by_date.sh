#!/bin/sh
# calls workouts with a particular version of matlab depending on the date today

N=4
dt=$(date +\%d)
dt_modN=$(($dt%N))

if [ $dt_modN == 0 ]
then
  matlabVer_1="R2017b"
  matlabVer_2="R2018a"
elif [ $dt_modN == 1 ]
then
  matlabVer_1="R2018b"
  matlabVer_2="R2019a"
elif [ $dt_modN == 2 ]
then
  matlabVer_1="R2019b"  # takes more than 20 hours
  # matlabVer_2="R2019b"
else
  matlabVer_1="R2020a"  # takes more than 20 hours
  # matlabVer_2="R2019b"
  #matlabVer_2="R2020b"
fi

echo "date" $dt_modN "matlab version" $matlabVer

sh /home/gail/GAIL_tests/PBS_jobs/gitpull.sh

# call the workouts
sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_workouts_report.sh $matlabVer_1

sh /home/gail/GAIL_tests/PBS_jobs/gitpull.sh

if [ -z "$matlabVer_2" ]
then
  echo "\$matlabVer_2 is empty"
else
  sh /home/gail/GAIL_tests/PBS_jobs/gailpbs_workouts_report.sh $matlabVer_2
fi
