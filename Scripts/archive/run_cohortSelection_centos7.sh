#$ -S /bin/sh   
#$ -M myuni@cumc.columbia.edu
#$ -m bea
#$ -N syn


PROJECT=$1
rfile=${PROJECT}/ClusteredCollapsing/Scripts/cohortSelection.R
project=${PROJECT}/ClusteredCollapsing/Results
sample_path="$PROJECT/cohort_in_dragen.txt"

# --dir ${dir}
source /nfs/goldstein/software/centos7/R-4.1.0_with_gcc_10-x86_64/R-4.1.0_with_gcc_10-x86_64-ENV.sh
R --slave --no-restore --file=${rfile} --args --case_list_path $sample_path 
