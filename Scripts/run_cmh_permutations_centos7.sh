#$ -S /bin/sh   
#$ -M myuni@cumc.columbia.edu
#$ -m bea
#$ -N syn

source /nfs/goldstein/software/centos7/python-3.9.7-x86_64_shared/python3.9.7-ENV.sh
PROJECT=$(cat ./Input/collapsing.yaml | shyaml get-value USER_VARIABLE.PROJECT)
# PROJECT=/nfs/projects/refractory_epilepsy/20220211_total/ClusteredCollapsing
rfile=${PROJECT}/ClusteredCollapsing/Scripts/cmh_permute_lclust_short_script.R
project=${PROJECT}/ClusteredCollapsing/Results
dir=""
res=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
# res="0_2"
minsample=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.minsample)
echo "Min Sample is $minsample"
# minsample=20
# model=dominantUltraRareEnsemble
model=URPTV_pext9_igm_af,URSyn_pext9_igm_af,FlexPTV_pext9 #URMisP9,URPTVMisP9,URPTVP9,domSynP9
echo $model
permstart=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.permstart)
echo "permstart is $permstart"
# permstart=1
permend=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.permend)
# permend=100
cores=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.cores)
# cores=12

# --dir ${dir}

/nfs/goldstein/software/centos7/R-4.0.0-x86_64/lib64/R/bin/R --slave --no-restore --file=${rfile} --args --project ${project}  --res ${res} --minsample ${minsample} --model ${model} --permstart ${permstart} --permend ${permend} --cores ${cores}
