# adapt this
PROJECT=$1
logdir=$PROJECT/ClusteredCollapsing

qsub -wd $logdir -q centos78.q -pe threaded 1 -V $PROJECT/ClusteredCollapsing/Scripts/run_cohortSelection_centos7.sh $1

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations_genderstrat_centos7.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations_genderstrat.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations_gender_centos7.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations_gender.sh
