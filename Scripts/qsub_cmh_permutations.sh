# adapt this
PROJECT=/nfs/projects/refractory_epilepsy/20220211_total/ClusteredCollapsing
logdir=$PROJECT/Results/Log

qsub -wd $logdir -q centos78.q -pe threaded 12 -V $PROJECT/Scripts/run_cmh_permutations_centos7.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations_genderstrat_centos7.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations_genderstrat.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations_gender_centos7.sh

# qsub -wd $logdir -pe threaded 24 -V $PROJECT/Scripts/run_cmh_permutations_gender.sh
