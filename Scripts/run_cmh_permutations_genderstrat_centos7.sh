#$ -S /bin/sh   
#$ -M myuni@cumc.columbia.edu
#$ -m bea
#$ -N RX

PROJECT=/pathtoproject/Collapsing1
rfile=${PROJECT}/Scripts/cmh_permute_lclust_short_genderstrat_script.R
project=${PROJECT}/Results
dir=All
res="0_3"
minsample=5
model=recessiveXwoMAPIDHP

permstart=1
permend=250
cores=24

/nfs/goldstein/software/centos7/R-4.0.0-x86_64/lib64/R/bin/R  --slave --no-restore --file=${rfile} --args --project ${project} --dir ${dir} --res ${res} --minsample ${minsample} --model ${model} --permstart ${permstart} --permend ${permend} --cores ${cores}
