#$ -S /bin/sh   
#$ -M myuni@cumc.columbia.edu
#$ -m bea
#$ -N RX

PROJECT=/pathtoproject/Collapsing1
rfile=${PROJECT}/Scripts/cmh_permute_lclust_short_gender_script.R
project=${PROJECT}/Results
dir=All
res="0_3"
minsample=5
gender="male"
# gender="female"
model=recessiveXwoMAPIDHP

permstart=1
permend=250
cores=24


/nfs/goldstein/software/R-3.4.3-x86_64/lib64/R/bin/R --slave --no-restore --file=$rfile --args --project ${project} --dir ${dir} --res ${res} --minsample ${minsample} --gender ${gender} --model ${model} --permstart ${permstart} --permend ${permend} --cores ${cores}
