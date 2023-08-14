#\!/bin/bash
# 2021-05-20, last updated 2023-04-14
# Run Collapsing Analysis 

time_stamp=$(date +%Y%m%d_%H%M%S)
################ Server specific ###################
source /nfs/goldstein/software/centos7/python-3.9.7-x86_64_shared/python3.9.7-ENV.sh
if [[ $1 != "runNHC" ]] && [[ $1 != "runNHCMixedAncestry" ]] 
then
  source /nfs/goldstein/software/centos7/R-4.1.0_with_gcc_10-x86_64/R-4.1.0_with_gcc_10-x86_64-ENV.sh
fi

PROJECT=$(pwd)
echo "The project folder is "
echo $PROJECT

atav="/nfs/goldstein/software/sh/atav.sh --email" # remove --email if you don't want to get them
atavtrunk="/nfs/goldstein/software/sh/atav_trunk.sh --email" # remove --email if you don't want to get them
Rscript_var="Rscript"
gene_bound_file_original_repo="/nfs/goldstein/software/atav_home/data/ccds/addjusted.CCDS.genes.index.r20.hg19.txt" #this is the file that was originally implemented in the repo
gene_bound_file_update_2022_1="/nfs/goldstein/software/atav_home/data/ccds/addjusted.CCDS.genes.index.r20.hg19.r15names.txt" #as of 7/18/2022, discussion that there had been previous attempts to use an updated gene boundary file with updated gene symbol names. Explanation here https://redmine.igm.cumc.columbia.edu/issues/3278
gene_bound_file_current="/nfs/goldstein/software/atav_home/data/ccds/addjusted.CCDS.genes.index.r20.hg19.ensembl87.txt" #currected. see https://redmine.igm.cumc.columbia.edu/issues/7377
regenie="/nfs/goldstein/software/centos7/regenie-3.1.3/regenie_v3.1.3.gz_x86_64_Centos7_mkl"
plink2="/usr/local/igm/non-atav-tools/plink2_20220814/plink2"

################ Paths to R scripts ###################
cohortSelection_var=$PROJECT/Scripts/cohortSelection.R
lclust_flash_var=$PROJECT/Scripts/lclust_Flash.R
lclust_flash_combo_cluster_var=$PROJECT/Scripts/lclust_Flash_combo_cluster.R
cmh_lclust_var=$PROJECT/Scripts/cmh_lclust.R
cmh_qq_var=$PROJECT/Scripts/cmh_qq.R
cmh_perm=$PROJECT/Scripts/cmh_permute_lclust_short_script.R
umap_clusters_fig=$PROJECT/Scripts/plot_umap_table.R
fp_forest_plot_var=$PROJECT/Scripts/fp_forest_plot.R
run_NHC_var=$PROJECT/Scripts/run_NHC.R
reg_covphe=$PROJECT/Scripts/regenieCovPheno.R
reg_annomaker=$PROJECT/Scripts/AnnoMaker.sh
reg_results=$PROJECT/Scripts/regenieResults.R
reg_results_SPA=$PROJECT/Scripts/regenieResults_SPA.R
reg_covphe_split=$PROJECT/Scripts/regenieCovPhenoSplit.R
reg_covphe_update=$PROJECT/Scripts/regenieCovPheno_update.R
power_collapsing_var=$PROJECT/Scripts/power_collapsing.R
unbiased_threshold_var=$PROJECT/Scripts/unbiased_threshold.R
parse_cov_detais_var=$PROJECT/Scripts/parseCoverageSummaryFileV2.pl
perform_coverage_pca_var=$PROJECT/Scripts/performCoveragePCA.R

if [[ -z $1 ]]
then
  echo "needs a parameter to determine which step to run"
  exit
fi

##################################################################################
# gives user feedback on whether bash argument was valid or not
##################################################################################
function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    LIST_WHITESPACES=`echo $LIST | tr "$DELIMITER" " "`
    for x in $LIST_WHITESPACES; do
        if [ "$x" = "$VALUE" ]; then
            return 0
        fi
    done
    return 1
}

ok_arg_list="cohortSelection initialCoverage CohortPCA clustering clusterCoverage clusterCoverageMale clusterCoverageFemale clusterCoverage clusterCoverageMale clusterCoverageFemale createMasterGeno models cmh qq umapTable sampleTracking RegListVCF clusterCoverageMixedAncestry clusteringMixedAncestry forestPlot_display pcaCluster pcaMixedAncestry summarizeRemovedSamples cmhPerm runNHC modelsMixedAncestry RegListVCFMixedAncestry RegPLINKMergeMixedAncestry_1 RegPLINKMergeMixedAncestry_2 RegCovPheMixedAncestry_1 RegCovPheMixedAncestry_2 RegCovPheMixedAncestry_split RegStep1MixedAncestryBT RegStep1MixedAncestryQT RegStep2ExWASBT RegStep2ExWASQT RegStep2VCtestsBT RegStep2VCtestsQT collapsingPower SampleFileCheck createMasterGenoMixedAncestry findThreshold TableWModels runNHCMixedAncestry RegStep2VCtestsBT_QQ create_digenic create_gene_set_model filter_loeuf RegStep2VCtests_QQManhattan RegStep2VCtests_QQManhattan_SPA RegStep2ExWASQT_SPA RegStep2ExWASBT_SPA RegStep2VCtestsBT_SPA RegStep2VCtestsQT_SPA combine_cluster_genotype geneAnno caseinModel RegListVCFMixedAncestry_genotype RegCovPheMixedAncestry_3 RegListVCFMixedAncestry_create_AAF RegStep2VCtestsBT_ch RegStep2VCtestsBT_max RegStep2VCtestsBT_sum"

if exists_in_list "$ok_arg_list" " " $1; then
  echo "$1 is a valid argument"
else
  echo "$1 is not a valid option. Valid options are $ok_arg_list"
  exit
fi

################ Read Me ###################
# Each step is run separately
# See documentation to determine location of code run (dev machines vs. qs1, future code will all run on qs1)

################ ###################
case_group=$(cat ./Input/collapsing.yaml | shyaml get-value USER_VARIABLE.case_group)
echo "case_group is $case_group"
################ run cohortSelection.R to get input file ###################

if [[ $1 = "cohortSelection" ]]
then

  if [[ $2 = "--help" ]]
  then
    # Learn about the input options for cohortSelection.R
      $Rscript_var $cohortSelection_var --help  
  else
      mkdir Results
      qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/"$time_stamp"_stdout_cohortSelection.log -e $PROJECT/Results/"$time_stamp"_stderr_cohortSelection.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "cohortSelection" $PROJECT
  fi
  exit
fi

###################################################################################################
###################################################################################################
#####                               Check that sample file is OK                              #####
###################################################################################################
###################################################################################################
SAMPLES="$PROJECT/Data/*.ped.txt"
if [[ $1 = "SampleFileCheck" ]]
then
  echo "Checking to make sure sample file is correct..."
  $atav --list-var --sample $SAMPLES --out $PROJECT/Results/SampleFileCheck --rs-number rs79585140
  exit
fi

###################################################################################################
###################################################################################################
#####                        Review sample file check: run on dev1-dev4                          #####
###################################################################################################
###################################################################################################
#review samples that were removed from the checked sample file
if [[ $1 = "summarizeRemovedSamples" ]]
then
	echo "Comparing samples present in SampleFileCheck..."
  
  #determine what samples were removed from ped.txt and checked sample file
  qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/"$time_stamp"_stdout_summarizeRemovedSamples.log -e $PROJECT/Results/"$time_stamp"_stderr_summarizeRemovedSamples.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "summarizeRemovedSamples" $PROJECT
  exit  
fi

###################################################################################################
###################################################################################################
#####                               COHORT PRUNING: run on qs1                                         #####
###################################################################################################
###################################################################################################
SAMPLES="$PROJECT/Results/SampleFileCheck/*_existing.sample.txt"
if [[ $1 = "initialCoverage" ]]
then
	echo "Creating Initial Coverage..."

	$atav --coverage-comparison --gene-boundary $gene_bound_file_current --min-coverage 10 --sample $SAMPLES --out $PROJECT/Results/Coverage
  exit
fi

##################################################################################
#####              KINSHIP + FlashPCA: run on qs1                           #####
##################################################################################
if [[ $1 = "CohortPCA" ]] 
then
  #this step is called "pruning" but this is a misnomer. Really, it just generates the PCAs which will later be used for clustering. Future directions should then prune within each cluster.  
	echo "Creating PCAs and pruning related individuals..."

  kinship_relatedness_threshold=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.kinship_relatedness_threshold)
  #Note: by default --flashpca will not perform PLINK nearest-neighbor pruning anymore use --flashpca-plink-pruning to trigger running plink pruning step if you want and add --flashpca-num-nearest-neighbor 200 
  $atav --ped-map --variant /nfs/goldstein/software/atav_home/data/variant/informative_snps.ld_pruned.37MB.txt --kinship --kinship-relatedness-threshold $kinship_relatedness_threshold --flashpca --min-covered-case-percentage 95 --min-covered-ctrl-percentage 95 --min-coverage 10 --sample-coverage-summary $PROJECT/Results/Coverage/*_sample.summary.csv --sample $SAMPLES --exclude-igm-gnomad-sample --out $PROJECT/Results/KinshipFlashPCA
  exit
fi

###############################################################################
#####                          Louvain Clustering: Run on qs1            ######
###############################################################################
if [[ $1 = "clustering" ]]
then

  if [[ $2 = "--help" ]]
  then
    # Learn about the input options for lclust_Flash.R
    $Rscript_var $lclust_flash_var --help 
  else
  	echo "Creating clusters..."
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/"$time_stamp"_stdout_clustering.log -e $PROJECT/Results/"$time_stamp"_stderr_clustering.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "clustering" $PROJECT
  fi
  exit
fi

# resolution="0_2" #this is defined by user
resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
echo "The resolution is $resolution"
sleep 2


##################################################################################
#####             Create a multi ancestry cluster: Run on Dev1-4            ######
##################################################################################

min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
echo "The min_sample is $min_sample"
max_ratio=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.max_ratio)
echo "The max ratio is $max_ratio"
if [[ $1 = "clusteringMixedAncestry" ]]
then
  if [[ $2 = "--help" ]]
  then
    # Learn about the input options for lclust_flash_combo_cluster.R
    $Rscript_var $lclust_flash_combo_cluster_var --help 
  else
  	echo "Create mixed ancestry cluster..."

    # run cluster creation
    # $Rscript_var $lclust_flash_combo_cluster_var --resolution_var $resolution --case_group $case_group --min_sample $min_sample --max_ratio $max_ratio
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/"$time_stamp"_stdout_clusteringMixedAncestry.log -e $PROJECT/Results/"$time_stamp"_stderr_clusteringMixedAncestry.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "clusteringMixedAncestry" $PROJECT
  fi
  exit
fi

##################################################################################
##### PCA by cluster: run on qs1. Not necessary as part of basic clustering but
##### provides cluster based PCAs for ease
##################################################################################

pca_start="$atav --ped-map --variant /nfs/goldstein/software/atav_home/data/variant/informative_snps.ld_pruned.37MB.txt --flashpca --flashpca-plink-pruning --flashpca-num-nearest-neighbor 200 "

min_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_cluster)
max_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.max_cluster)
echo "The min_cluster is $min_cluster"
echo "The max_cluster is $max_cluster"
sleep 2

if [[ $1 = "pcaCluster" ]]
then
	echo "Running PCA per cluster..."

  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do
 
  samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
  out="$PROJECT/Results/KinshipFlashPCA/PCAFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07/"
  
  # echo $i
  
  pca_end=" --sample $samples --out $out"
  pca_clust_exec_str="$pca_start$pca_end"
  echo $pca_clust_exec_str
  eval $pca_clust_exec_str
  sleep 2
  done
  exit
fi

if [[ $1 = "pcaMixedAncestry" ]]
then
	echo "Running pca for mixed ancestry_cluster..."

  mixed_ancestry_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.mixed_ancestry_cluster)
  echo "The mixed ancestry cluster is $mixed_ancestry_cluster"
	echo $max_ratio

  samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/*flashPCA_lclustering_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_sample.txt"
  out="$PROJECT/Results/KinshipFlashPCA/PCAFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/"

  pca_end=" --sample $samples --out $out"
  pca_clust_exec_str="$pca_start$pca_end"
  echo $pca_clust_exec_str
  eval $pca_clust_exec_str

  exit
fi

##################################################################################
#####                       COVERAGE2: Run on qs1                           ######
##################################################################################
cov_clust_exec_str_start="$atav --site-coverage-comparison --gene-boundaries $gene_bound_file_current --site-max-percent-cov-difference 0.07 --min-coverage 10"

excludeOutliers=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.excludeOutliers)
# echo $excludeOutliers

if [[ $1 = "clusterCoverage" ]]
then
	echo "Running coverage per cluster..."
  usecoveragePCA=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.useCoveragePCA)
  runcoveragePCA=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.runCoveragePCA)

  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do 

    if [[ $usecoveragePCA = "True" ]]
    then
      echo "usecoveragePCA is TRUE"
      samples="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_CovPCA/*sd3_filtered.txt"
    elif [[ $excludeOutliers = "True" ]] 
    then
      echo "excludeOutliers is TRUE"
      samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"

      lin_num=$(cat $samples | wc -l)
      echo "there are $lin_num samples"
      if [[ $lin_num -gt 200 ]] #pruning only works with at least 200 samples
      then
        echo "Therefore looking for pruned sample file"
        samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/PCAFlashLClust_res_"${resolution}"_cluster_"${i}"_07/*lashpca_pruned_sample_file.txt"
      else
        echo "Therefore NOT looking for pruned sample file"
        samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
      fi    
    elif [[ $excludeOutliers = "False" ]]
    then
      echo "excludeOutliers is FALSE"
      samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
    else
      echo "excludeOutliers is $excludeOutliers It must be TRUE or FALSE"
      exit
    fi
    out="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07/"

    if [[ $runcoveragePCA = "True" ]] && [[ $usecoveragePCA = "False" ]]
    then
      if [[ -d "$PROJECT/Results/KinshipFlashPCA"${dir}"/PCAFlashLClust_res_"${resolution}"_cluster_"${i}"_07/" ]]
      then
        echo "coveragePCA is TRUE"

        cov_clust_exec_str_start="$atav --site-coverage-comparison --gene-boundaries $gene_bound_file_current --site-max-percent-cov-difference 0.07 --min-coverage 10 --include-coverage-detail"
        out="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_CovPCA/"
      else
        echo "Must first run pcaCluster before you can runcoveragePCA even if you will not exclude outliers"
        exit
      fi
    fi

    cov_clust_exec_str_end=" --sample $samples --out $out"
    cov_clust_exec_str="$cov_clust_exec_str_start$cov_clust_exec_str_end"
    echo $cov_clust_exec_str
    eval $cov_clust_exec_str
    sleep 5
    done
    exit
fi

if [[ $1 = "clusterCoverageMale" ]]
then
	echo "Running coverage per cluster for Males..."

  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do
  
  samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample_male.txt"
  out="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_male/"
  
  cov_clust_exec_str_end=" --sample $samples --out $out"
  cov_clust_exec_str="$cov_clust_exec_str_start$cov_clust_exec_str_end"
  echo $cov_clust_exec_str
  eval $cov_clust_exec_str
  sleep 2
  
  done
  exit
fi

if [[ $1 = "clusterCoverageFemale" ]]
then
	echo "Running coverage per cluster for Females..."

  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do
  
  samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample_female.txt"
  out="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_female/"
  
  cov_clust_exec_str_end=" --sample $samples --out $out"
  cov_clust_exec_str="$cov_clust_exec_str_start$cov_clust_exec_str_end"
  echo $cov_clust_exec_str
  eval $cov_clust_exec_str
  sleep 2
  done
  exit
fi

if [[ $1 = "clusterCoverageMixedAncestry" ]]
then
	echo "Running coverage for mixed ancestry_cluster..."

  mixed_ancestry_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.mixed_ancestry_cluster)
  echo "The resolution is $mixed_ancestry_cluster"

  if [[ $excludeOutliers = "True" ]] 
  then
    echo "excludeOutliers is TRUE"
    samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/*flashPCA_lclustering_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_sample.txt"

    lin_num=$(cat $samples | wc -l)
    echo "there are $lin_num samples"
    if [[ $lin_num -gt 200 ]] #pruning only works with at least 200 samples
    then
      echo "Therefore looking for pruned sample file"
      samples="$PROJECT/Results/KinshipFlashPCA/PCAFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio_label}"_07_mixedCluster/*lashpca_pruned_sample_file.txt"
    else
      echo "Therefore NOT looking for pruned sample file"
      samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/*flashPCA_lclustering_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_sample.txt"
    fi    
  elif [[ $excludeOutliers = "False" ]]
  then
    echo "excludeOutliers is FALSE"
    samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/*flashPCA_lclustering_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_sample.txt"
  else
    echo "excludeOutliers is $excludeOutliers It must be TRUE or FALSE"
    exit
  fi

  out="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/"

  cov_clust_exec_str_end=" --sample $samples --out $out"
  cov_clust_exec_str="$cov_clust_exec_str_start$cov_clust_exec_str_end"
  # echo $out

  echo $cov_clust_exec_str
  eval $cov_clust_exec_str
  exit
fi

###################################################################################################
###################################################################################################
#####                                   COLLAPSING: Run on qs1                                #####
###################################################################################################
###################################################################################################

CODING_SPLICE="HIGH:exon_loss_variant,HIGH:frameshift_variant,HIGH:rare_amino_acid_variant,HIGH:stop_gained,HIGH:start_lost,HIGH:stop_lost,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant,HIGH:gene_fusion,HIGH:bidirectional_gene_fusion,MODERATE:3_prime_UTR_truncation+exon_loss_variant,MODERATE:5_prime_UTR_truncation+exon_loss_variant,MODERATE:coding_sequence_variant,MODERATE:disruptive_inframe_deletion,MODERATE:disruptive_inframe_insertion,MODERATE:conservative_inframe_deletion,MODERATE:conservative_inframe_insertion,MODERATE:missense_variant+splice_region_variant,MODERATE:missense_variant,MODERATE:splice_region_variant,LOW:5_prime_UTR_premature_start_codon_gain_variant,LOW:initiator_codon_variant,LOW:initiator_codon_variant+non_canonical_start_codon,LOW:splice_region_variant+synonymous_variant,LOW:splice_region_variant,LOW:start_retained,LOW:stop_retained_variant,LOW:synonymous_variant"

FUNCTIONAL_EFFECTS="HIGH:exon_loss_variant,HIGH:frameshift_variant,HIGH:rare_amino_acid_variant,HIGH:stop_gained,HIGH:start_lost,HIGH:stop_lost,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant,HIGH:gene_fusion,HIGH:bidirectional_gene_fusion,MODERATE:3_prime_UTR_truncation+exon_loss_variant,MODERATE:5_prime_UTR_truncation+exon_loss_variant,MODERATE:coding_sequence_variant,MODERATE:disruptive_inframe_deletion,MODERATE:disruptive_inframe_insertion,MODERATE:conservative_inframe_deletion,MODERATE:conservative_inframe_insertion,MODERATE:missense_variant+splice_region_variant,MODERATE:missense_variant,LOW:5_prime_UTR_premature_start_codon_gain_variant,LOW:initiator_codon_variant,LOW:initiator_codon_variant+non_canonical_start_codon"

MOD_HIGH_EFFECTS="HIGH:exon_loss_variant,HIGH:frameshift_variant,HIGH:rare_amino_acid_variant,HIGH:stop_gained,HIGH:start_lost,HIGH:stop_lost,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant,HIGH:gene_fusion,HIGH:bidirectional_gene_fusion,MODERATE:3_prime_UTR_truncation+exon_loss_variant,MODERATE:5_prime_UTR_truncation+exon_loss_variant,MODERATE:coding_sequence_variant,MODERATE:disruptive_inframe_deletion,MODERATE:disruptive_inframe_insertion,MODERATE:conservative_inframe_deletion,MODERATE:conservative_inframe_insertion,MODERATE:missense_variant+splice_region_variant,MODERATE:missense_variant,MODERATE:splice_region_variant"

FUNCTIONAL_EFFECTS_HQ="HIGH:frameshift_variant,HIGH:stop_gained,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant,MODERATE:missense_variant+splice_region_variant,MODERATE:missense_variant"

LOF_EFFECTS="HIGH:exon_loss_variant,HIGH:frameshift_variant,HIGH:rare_amino_acid_variant,HIGH:stop_gained,HIGH:stop_lost,HIGH:start_lost,HIGH:gene_fusion,HIGH:bidirectional_gene_fusion,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant"

LOF_EFFECTS_HQ="HIGH:frameshift_variant,HIGH:stop_gained,HIGH:splice_acceptor_variant,HIGH:splice_donor_variant"

SYN_EFFECTS="LOW:start_retained,LOW:stop_retained_variant,LOW:synonymous_variant"

MISSENSE_ONLY="MODERATE:missense_variant+splice_region_variant,MODERATE:missense_variant"

INCLUDES="--include-gnomad-genome --include-gnomad-exome --include-mtr --include-gerp --include-rvis --include-sub-rvis --include-limbr --include-revel --include-trap --include-known-var --include-primate-ai --include-ccr --include-loftee --include-gnomad-gene-metrics --include-pext --include-mpc --flag-repeat-region --include-syn-rvis  --include-top-med --include-igm-af" # --include-discovehr --include-genome-asia --include-iranome --include-gme

QC="--filter pass,likely,intermediate --ccds-only --min-coverage 10 --include-qc-missing --qd 5 --qual 50 --mq 40 --gq 20 --snv-fs 60 --indel-fs 200 --snv-sor 3 --indel-sor 10 --rprs -3 --mqrs -10 --het-percent-alt-read 0.3-1 --gnomad-genome-rf-tp-probability-snv 0.01 --gnomad-genome-rf-tp-probability-indel 0.02 --gnomad-exome-rf-tp-probability-snv 0.01 --gnomad-exome-rf-tp-probability-indel 0.02"

GNOMAD_NONNEURO_POP="non_neuro_afr,non_neuro_amr,non_neuro_asj,non_neuro_eas,non_neuro_sas,non_neuro_fin,non_neuro_nfe"

GNOMAD_POP="afr,amr,asj,eas,sas,fin,nfe"

GNOMAD_NONNEURO_POP_GENOME="non_neuro_afr,non_neuro_amr,non_neuro_asj,non_neuro_eas,non_neuro_fin,non_neuro_nfe"
GNOMAD_POP_GENOME="afr,amr,asj,eas,fin,nfe"

INCLUDES_ATAVTRUNK="" 

AUTOSOMES="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"

#################################################################################################


#################################################################################################
############# Run Master genotypes to be used for collapsing lite later: Run on qs1
#################################################################################################
if [[ $1 = "createMasterGeno" ]]
then
  echo "Creating master genotype.csv ..."
  
  
  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do 
  outputFolder="$PROJECT/Results/Collapsing"${dir}"/LClust_res_"${resolution}"_cluster_"${i}_${dir2}"FlashColl_07"
  echo $samples
  
  
  # gundi version: max-qc-fail-sample missing, lenient effect filter, loose AF filter
  # $atav --collapsing-dom --mann-whitney-test --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES $QC --effect $CODING_SPLICE --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop global --gnomad-exome-maf 0.01 --exac-pop global --exac-maf 0.01 --loo-maf 0.01 --sample $samples --out $outputFolder/dominantNoneMAF/
  # sleep 2

    if [[ -z $2 ]] #using site-based coverage pruning. will use whatever sample file was used to generate the coverage.
    then
      # Josh does same version but with list-var geno to save a little sapce
      samples="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}"_07/*existing.sample.txt"
      geneBoundaries="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07/*_site.clean.txt"
      coverage2="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07/*_coverage.summary.clean.csv"
      $atav --list-var-geno --gene-boundaries $geneBoundaries $INCLUDES $QC --effect $CODING_SPLICE --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop global --gnomad-exome-maf 0.01 --loo-maf 0.01 --gzip --sample $samples --out $outputFolder/dominantNoneMAF/
      sleep 2
    elif [[ $2 = "skipCoverage" ]]  #skipping site based coverage pruning.
    then
      if [[ $excludeOutliers = "False" ]] 
      then
        echo "excludeOutliers is FALSE"
        samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
      elif [[ $excludeOutliers = "True" ]] 
      then
        echo "excludeOutliers is TRUE"
        samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"

        lin_num=$(cat $samples | wc -l)
        echo "there are $lin_num samples"
        if [[ $lin_num -gt 200 ]] #pruning only works with at least 200 samples
        then
          echo "Therefore looking for pruned sample file"
          samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/PCAFlashLClust_res_"${resolution}"_cluster_"${i}"_07/*lashpca_pruned_sample_file.txt"
        else
          echo "Therefore NOT looking for pruned sample file"
          samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
        fi
      else
        echo "excludeOutliers is $excludeOutliers It must be TRUE or FALSE"
        exit
      fi      
      $atav --list-var-geno --gene-boundaries $gene_bound_file_current $INCLUDES $QC --site-max-percent-cov-difference 0.07 --effect $CODING_SPLICE --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop global --gnomad-exome-maf 0.01 --loo-maf 0.01 --gzip --sample $samples --out $outputFolder/dominantNoneMAF/
        sleep 2
    else
      echo "Second argument not supported"
      exit
    fi
  done
  exit
fi

if [[ $1 = "createMasterGenoMixedAncestry" ]]
then
  echo "Creating master genotype file for mixed ancestry_cluster..."
  sleep 2

  mixed_ancestry_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.mixed_ancestry_cluster)
  echo "The mixed ancestry cluster is $mixed_ancestry_cluster"
  
  outputFolder="$PROJECT/Results/Collapsing"${dir}"/LClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_mixedCluster_FlashColl_07/"

  if [[ -z $2 ]] #using site-based coverage pruning
  then
    samples=$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/*_existing.sample.txt
    geneBoundaries="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/*_site.clean.txt"
    $atav --list-var-geno --gene-boundaries $geneBoundaries $INCLUDES $QC --effect $CODING_SPLICE --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop global --gnomad-exome-maf 0.01  --loo-maf 0.01 --sample $samples --out $outputFolder/dominantNoneMAF/
  elif [[ $2 = "skipCoverage" ]] #skipping site based coverage pruning.
  then
    echo "Skipping coverage. At the moment, this step uses the sample file from the clustering step. It does not use the sample file from any PCA pruning steps run separately on the mixed ancestry cluster."
    sleep 10
    samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/*flashPCA_lclustering_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio_label}"_sample.txt"
    $atav --list-var-geno --gene-boundaries $gene_bound_file_current $INCLUDES $QC --site-max-percent-cov-difference 0.07 --effect $CODING_SPLICE --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop global --gnomad-exome-maf 0.01 --loo-maf 0.01 --sample $samples --out $outputFolder/dominantNoneMAF/
    else
      echo "Second argument not supported"
  fi

  exit

  cov_clust_exec_str_end=" --sample $samples --out $out"
  cov_clust_exec_str="$cov_clust_exec_str_start$cov_clust_exec_str_end"
  # echo $out

  echo $cov_clust_exec_str
  eval $cov_clust_exec_str
  exit
fi

#################################################################################################
############# Create Models 
#################################################################################################
# if models or mixed ancestry models, get models

if [[ $1 = "models" ]] || [[ $1 = "modelsMixedAncestry" ]] #if just models, execute on min to mx
then
  echo "Creating models..."

  model_num=$(cat ./Input/collapsing.yaml | shyaml get-value MODEL_VARIABLE.model_num)
  echo "The Model numbers to run are $model_num"

  list_models=()
  
  #list all models to run in new array
  for x in ${model_num[@]}
  do
    #if x contains -, loop through start and end num
    if [[ $x == *"-"* ]]
    then
      Start=$(echo ${x/-/" "} | awk '{print $1}')    
      End=$(echo ${x/-/" "} | awk '{print $2}')      
      for (( c=$Start; c<=$End; c++ ))    
      do
        list_models+=( $c )
      done
    else
      list_models+=( $x )
    fi
  done 
fi
if [[ $1 = "models" ]] 
then

  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do 
  
    outputFolder="$PROJECT/Results/Collapsing"${dir}"/LClust_res_"${resolution}"_cluster_"${i}_${dir2}"FlashColl_07"
    genos=$outputFolder"/dominantNoneMAF/*genotypes.csv.gz"
    samples=$outputFolder"/dominantNoneMAF/*existing.sample.txt"
    coverage2="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07/*_coverage.summary.clean.csv"
    geneBoundaries="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07/*_site.clean.txt"
    
    for x in ${list_models[@]}
    do
      model=$(cat ./Input/collapsing.yaml | shyaml get-value MODEL_VARIABLE.model_${x})
      echo $model
      eval $model
      sleep 2
    done 

  done
  exit
fi

if [[ $1 = "modelsMixedAncestry" ]] #execute models on mixed ancestry model
then
  echo "Creating mixed ancesty models..."
  mixed_ancestry_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.mixed_ancestry_cluster)
  echo "The mixed ancestry cluster is $mixed_ancestry_cluster"

  outputFolder="$PROJECT/Results/Collapsing"${dir}"/LClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_mixedCluster_FlashColl_07/"

  genos=$outputFolder"/dominantNoneMAF/*genotypes.csv"
  samples=$outputFolder"/dominantNoneMAF/*existing.sample.txt"
  echo $genos

  for x in ${list_models[@]}
  do
    model=$(cat ./Input/collapsing.yaml | shyaml get-value MODEL_VARIABLE.model_${x})
    eval $model
    sleep 2
  done 
  exit
fi

if [[ $1 = "not_active" ]]
then

  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do 
  samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
  geneBoundaries="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07/*_site.clean.txt"
  coverage2="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07/*_coverage.summary.clean.csv"
  outputFolder="$PROJECT/Results/Collapsing"${dir}"/LClust_res_"${resolution}"_cluster_"${i}_${dir2}"FlashColl_07"
  # echo $samples
  
  
  # sleep 2
  done
  
  
  
  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do 
  
  samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample_male.txt"
  geneBoundaries="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_male/*_site.clean.txt"
  coverage2="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_male/*_coverage.summary.clean.csv"
  outputFolder="$PROJECT/Results/Collapsing"${dir}"/LClust_res_"${resolution}"_cluster_"${i}_${dir2}"FlashColl_07_male"
  # echo $samples
  
  # $atav --collapsing-comp-het --exclude-comp-het-pid-variant --exclude-comp-het-hp-variant --mann-whitney-test --region "X" --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $FUNCTIONAL_EFFECTS --exclude-multiallelic-variant-2 $QC  --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop afr,amr,asj,eas,sas,fin,nfe --gnomad-exome-maf 0.01 --exac-pop afr,amr,nfe,fin,eas,sas --exac-maf 0.01 --loo-maf 0.01 --max-qc-fail-sample 2 --sample $samples --out $outputFolder/recessiveXwoMAPIDHP/
  # sleep 2
  
  # $atav --collapsing-comp-het --exclude-comp-het-pid-variant --exclude-comp-het-hp-variant --mann-whitney-test --region "X" --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $FUNCTIONAL_EFFECTS --polyphen probably,possibly,unknown --min-primate-ai 0.5 --min-revel-score 0.25 --ensemble-missense --exclude-multiallelic-variant-2 $QC  --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop afr,amr,asj,eas,sas,fin,nfe --gnomad-exome-maf 0.01 --exac-pop afr,amr,nfe,fin,eas,sas --exac-maf 0.01 --loo-maf 0.01 --max-qc-fail-sample 2 --sample $samples --out $outputFolder/recessiveXLEwoMAPIDHP/
  # sleep 2
  
  samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample_female.txt"
  geneBoundaries="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_female/*_site.clean.txt"
  coverage2="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_female/*_coverage.summary.clean.csv"
  outputFolder="$PROJECT/Results/Collapsing"${dir}"/LClust_res_"${resolution}"_cluster_"${i}_${dir2}"FlashColl_07_female"
  # echo $samples
  
  # $atav --collapsing-comp-het --exclude-comp-het-pid-variant --exclude-comp-het-hp-variant --mann-whitney-test --region "X" --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $FUNCTIONAL_EFFECTS --exclude-multiallelic-variant-2 $QC  --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop afr,amr,asj,eas,sas,fin,nfe --gnomad-exome-maf 0.01 --exac-pop afr,amr,nfe,fin,eas,sas --exac-maf 0.01 --loo-maf 0.01 --max-qc-fail-sample 2 --sample $samples --out $outputFolder/recessiveXwoMAPIDHP/
  # sleep 2
  
  # $atav --collapsing-comp-het --exclude-comp-het-pid-variant --exclude-comp-het-hp-variant --mann-whitney-test --region "X" --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $FUNCTIONAL_EFFECTS --polyphen probably,possibly,unknown --min-primate-ai 0.5 --min-revel-score 0.25 --ensemble-missense --exclude-multiallelic-variant-2 $QC  --gnomad-genome-pop global --gnomad-genome-maf 0.01 --gnomad-exome-pop afr,amr,asj,eas,sas,fin,nfe --gnomad-exome-maf 0.01 --exac-pop afr,amr,nfe,fin,eas,sas --exac-maf 0.01 --loo-maf 0.01 --max-qc-fail-sample 2 --sample $samples --out $outputFolder/recessiveXLEwoMAPIDHP/
  # sleep 2
  
  
  done
  
  
  
  # ### Dominant Syn
  # $atav --collapsing-dom --mann-whitney-test --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $SYN_EFFECTS $QC --min-exac-vqslod-snv 5000 --min-exac-vqslod-indel 5000  --gnomad-genome-pop global --gnomad-genome-maf 0 --gnomad-exome-pop global --gnomad-exome-maf 0 --exac-pop global --exac-maf 0 --loo-maf 0.0005 --max-qc-fail-sample 0 --sample $samples --out $outputFolder/dominantSynonymous/
  # sleep 2
  
  # ### Dominant Ultra-rare
  # $atav --collapsing-dom --mann-whitney-test --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $FUNCTIONAL_EFFECTS --polyphen probably $QC --min-exac-vqslod-snv 5000 --min-exac-vqslod-indel 5000  --gnomad-genome-pop global --gnomad-genome-maf 0 --gnomad-exome-pop global --gnomad-exome-maf 0 --exac-pop global --exac-maf 0 --loo-maf 0.0005 --max-qc-fail-sample 0 --sample $samples --out $outputFolder/dominantUltraRare/
  # sleep 2
  
  # ### Dominant PTV
  # $atav --collapsing-dom --mann-whitney-test --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $LOF_EFFECTS $QC  --gnomad-genome-pop global --gnomad-genome-maf 0.001 --gnomad-exome-pop afr,amr,asj,eas,sas,fin,nfe --gnomad-exome-maf 0.001 --exac-pop afr,amr,nfe,fin,eas,sas --exac-maf 0.001 --loo-maf 0.001 --max-qc-fail-sample 2 --sample $samples --out $outputFolder/dominantPTV/
  # sleep 2
  
  # ### Dominant Flexible Polyphen
  # $atav --collapsing-dom --mann-whitney-test --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $FUNCTIONAL_EFFECTS --polyphen probably $QC  --gnomad-genome-pop global --gnomad-genome-maf 0.001 --gnomad-exome-pop afr,amr,asj,eas,sas,fin,nfe --gnomad-exome-maf 0.001 --exac-pop afr,amr,nfe,fin,eas,sas --exac-maf 0.001 --loo-maf 0.001 --max-qc-fail-sample 2 --sample $samples --out $outputFolder/dominantFlexiblePolyphenDamaging/
  # sleep 2
  
  # ### Dominant Flexible #2
  # $atav --collapsing-dom --mann-whitney-test --gene-boundaries $geneBoundaries --read-coverage-summary $coverage2 $INCLUDES --effect $FUNCTIONAL_EFFECTS $QC  --gnomad-genome-pop global --gnomad-genome-maf 0.001 --gnomad-exome-pop afr,amr,asj,eas,sas,fin,nfe --gnomad-exome-maf 0.001 --exac-pop afr,amr,nfe,fin,eas,sas --exac-maf 0.001 --loo-maf 0.001 --max-qc-fail-sample 2 --sample $samples --out $outputFolder/dominantFlexible2/
  # sleep 2
fi

##################################################################################
#####                                   CMH                                 ######
##################################################################################

#########################################################################
#####                          run cmh_lclust.R  run on qs1        ######
#########################################################################
if [[ $1 = "cmh" ]]
then
  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $cmh_lclust_var --help 
  else
    echo "Creating cmh..."

    # run script
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/"$time_stamp"_stdout_cmh.log -e $PROJECT/Results/"$time_stamp"_stderr_cmh.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "cmh" $PROJECT
  fi 
  exit
fi

##################################################################################
#####                          run on qs1            ######
##################################################################################
if [[ $1 = "cmhPerm" ]]
then
  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $cmh_perm --help 
  else
    echo "Running cmh perm..."
    
    cores=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.cores)
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/"$time_stamp"_stdout_cmhPerm.log -e $PROJECT/Results/"$time_stamp"_stderr_cmhPerm.log  -pe smp $cores -V Scripts/run_collapsing_lclust.sh "cmhPerm" $PROJECT

  fi
  exit
fi


#####################################################################
#####                          run cmh_qq.R on qs1             ######
#####################################################################
if [[ $1 = "qq" ]]
then
  echo "Creating qq plots ..."
  mkdir $PROJECT/Results/CMH
  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $cmh_qq_var --help
  else
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/CMH/"$time_stamp"_stdout_qq.log -e $PROJECT/Results/CMH/"$time_stamp"_stderr_qq.log  -pe smp 1 -V Scripts/run_collapsing_lclust.sh "qq" $PROJECT
  fi
  exit
fi


##################################################################################
#####                          create gene_set             ######
##################################################################################
if [[ $1 = "create_gene_set_model" ]]
then

  create_gene_set=$PROJECT/Scripts/fp_create_gene_set.R
  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $create_gene_set --help
  else
    for (( i=$min_cluster; i<=$max_cluster; i++ ))
    do 
      qsub -wd $PROJECT -o $PROJECT/Results/Collapsing/"$time_stamp"_stdout_create_gene_set_model_cluster_"$i".log -e $PROJECT/Results/Collapsing/"$time_stamp"_stderr_create_gene_set_model_cluster_"$i".log -q centos78.q -pe smp 1 -V Scripts/run_collapsing_lclust.sh "create_gene_set_model" $PROJECT $i $gene_bound_file_current
      sleep 2
    done 
    sleep 90

    for (( i=$min_cluster; i<=$max_cluster; i++ ))
    do 
      

      resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
      gene_set_path=$(cat ./Input/collapsing.yaml | shyaml get-value GENE_SET.path)
      model=$(cat ./Input/collapsing.yaml | shyaml get-value GENE_SET.model)
      gene_set_name=$(cat ./Input/collapsing.yaml | shyaml get-value GENE_SET.name)

      model_folder="$PROJECT/Results/Collapsing/LClust_res_"${resolution}"_cluster_"${i}"_FlashColl_07"/$model
      echo "model folder"
      echo $model_folder
      outputFolder=$model_folder"_"$gene_set_name
      echo "outputfolder"
      echo $outputFolder
      genos="$outputFolder"/*genotypes1.csv""
      echo "genos"
      echo $genos
      samples="$model_folder"/*existing.sample.txt""
      echo "samples"
      echo $samples

      $atav --collapsing-lite --mann-whitney-test --genotype $genos --sample $samples --out $outputFolder
      sleep 2
    done 
  fi
  exit
fi

##################################################################################
#####                          filter loeuf             ######
##################################################################################
if [[ $1 = "filter_loeuf" ]]
then

  filter_loeuf=$PROJECT/Scripts/fp_filter_loeuf.R
  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $filter_loeuf --help
  else
    for (( i=$min_cluster; i<=$max_cluster; i++ ))
    do 
      resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
      model=$(cat ./Input/collapsing.yaml | shyaml get-value LOEUF_FILTER.model)
      loeuf_disp=$(cat ./Input/collapsing.yaml | shyaml get-value LOEUF_FILTER.loeuf | sed 's/\./_/')
      
      qsub -wd $PROJECT -o $PROJECT/Results/Collapsing/stdout.log -e $PROJECT/Results/Collapsing/stderr.log -q centos78.q -pe smp 1 -V Scripts/run_collapsing_lclust.sh "filter_loeuf" $PROJECT $i
      sleep 60
      
      model_folder="$PROJECT/Results/Collapsing/LClust_res_"${resolution}"_cluster_"${i}"_FlashColl_07"/$model
      echo "model folder"
      echo $model_folder
      outputFolder=$model_folder"_loeuf_"$loeuf_disp
      echo "outputfolder"
      echo $outputFolder
      genos="$outputFolder"/*genotypes1.csv""
      echo "genos"
      echo $genos
      samples="$model_folder"/*existing.sample.txt""
      echo "samples"
      echo $samples

      $atav --collapsing-lite --mann-whitney-test --genotype $genos --sample $samples --out $outputFolder
    done 
  fi
  exit
fi
##################################################################################
#####                          create digenic             ######
##################################################################################
if [[ $1 = "create_digenic" ]]
then

  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    create_digenic_var=$PROJECT/Scripts/digenic_create_model.R
    $Rscript_var $create_digenic_var --help
  else
    echo "Creating digenic plots ..."
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/stdout.log -e $PROJECT/Results/stderr.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "create_digenic" $PROJECT $2
  fi
  exit
fi

##################################################################################
#####                          combine cluster genotype files             ######
##################################################################################
if [[ $1 = "combine_cluster_genotype" ]]
then

  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    review_model=$PROJECT/Scripts/review_genotypes_from_model.R
    $Rscript_var $review_model --help
  else
    echo "combine genotype files ..."
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/"$time_stamp"_stdout_combine_cluster_genotype.log -e $PROJECT/Results/"$time_stamp"_stderr_combine_cluster_genotype.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "combine_cluster_genotype" $PROJECT
  fi
  exit
fi

##################################################################################
#####                          run plot_umap_table.R on dev1 to dev4        ######
##################################################################################
if [[ $1 = "umapTable" ]]
then
  echo "Creating umap with clusters and tables plots ..."
  
  # Learn about the input options
  # $Rscript_var $umap_clusters_fig --help
  
  # run script
  # $Rscript_var $umap_clusters_fig --resolution_var $resolution --min_sample 20 --case_group $case_group
fi

##########################################
#####             run on qs1        ######
##########################################
if [[ $1 = "sampleTracking" ]]
then
  sample_counts_table=$PROJECT/Scripts/table_sample_counts.R
  if [[ $2 = "--help" ]]
  then

    # Learn about the input options
    $Rscript_var $sample_counts_table --help
  else
    echo "Creating sample tracking log..."
    
    # run script
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/"$time_stamp"_sampleTracking_stdout.log -e $PROJECT/Results/"$time_stamp"_sampleTracking_stderr.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "sampleTracking" $PROJECT
  fi  
  exit
fi


##################################################################################
#####                          run tableWModels.R on dev machines    ######
##################################################################################
if [[ $1 = "TableWModels" ]]
then
  echo "Creating table with models..."
  if [[ $2 = "--help" ]]
  then
    table_w_models=$PROJECT/Scripts/TableWModels.R
    # Learn about the input options
    $Rscript_var $table_w_models --help
  else
    mkdir $PROJECT/Results/TableWModels
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/TableWModels/"$time_stamp"_TableWModels_stdout.log -e $PROJECT/Results/TableWModels/"$time_stamp"_TableWModels_stderr.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "TableWModels" $PROJECT
  fi  
  exit
fi


###################################################################################
# Creates forest plot. Run on qs1
##################################################################################
if [[ $1 = "forestPlot_display" ]]
then
  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $fp_forest_plot_var --help
  else
    echo "Creating forest plot..."
    mkdir $PROJECT/Results/fp
    cores=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.cores)
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/fp/"$time_stamp"_stdout_forestPlot_display.log -e $PROJECT/Results/fp/"$time_stamp"_stderr_forestPlot_display.log -pe smp 2 -V Scripts/run_collapsing_lclust.sh "forestPlot_display" $PROJECT
  fi
  exit
fi

##################################################################################
# Runs power/odds ratio calc Run on dev1:dev4
##################################################################################
if [[ $1 = "collapsingPower" ]]
then
  calc=$(cat ./Input/collapsing.yaml | shyaml get-value POWER_COLLAPSING.calc)
  alpha=$(cat ./Input/collapsing.yaml | shyaml get-value POWER_COLLAPSING.alpha)
  odds_ratio=$(cat ./Input/collapsing.yaml | shyaml get-value POWER_COLLAPSING.odds_ratio)

  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $power_collapsing_var --help
  else
    echo "Running collapsing power..."

    # run script
    $Rscript_var $power_collapsing_var --calc $calc --odds_ratio $odds_ratio --alpha $alpha --case_group $case_group
  fi
  exit
fi

##################################################################################
# Runs coveragepca Run on dev1:dev4
##################################################################################
if [[ $1 = "pcaCoverage" ]]
then
  source /nfs/goldstein/software/centos7/perl-5.32-x86_64/perl-5.32-ENV.sh

  echo "Running coverage PCA.."

  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    perl $parse_cov_detais_var --help

    $Rscript_var $perform_coverage_pca_var --help
  else
    for (( i=$min_cluster; i<=$max_cluster; i++ ))
    do
      echo "Runnning cluster ${i}"

      if [[ $excludeOutliers = "True" ]] 
      then
        echo "excludeOutliers is TRUE"
        samples="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"

        lin_num=$(cat $samples | wc -l)
        echo "there are $lin_num samples"
        if [[ $lin_num -gt 200 ]] #pruning only works with at least 200 samples
        then
          echo "Therefore looking for pruned sample file"
          sample_file="$PROJECT/Results/KinshipFlashPCA"${dir}"/PCAFlashLClust_res_"${resolution}"_cluster_"${i}"_07/*lashpca_pruned_sample_file.txt"
        else
          echo "Therefore NOT looking for pruned sample file"
          sample_file="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
        fi
      else
        echo "excludeOutliers is FALSE"

        sample_file="$PROJECT/Results/KinshipFlashPCA"${dir}"/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
      fi

      cov_details="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_CovPCA/*_coverage.detail.csv"
      perl $parse_cov_detais_var --cov_details $cov_details --sample_file $sample_file --cluster ${i}

      pca_file="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}_${dir2}"07_CovPCA/*.READYforPCA.txt"
      $Rscript_var $perform_coverage_pca_var --ready_for_pca_file $pca_file --cluster ${i}
    done
  fi
  exit
fi

#################################################################
#####              Gene Annotation                          #####
#################################################################
if [[ $1 = "geneAnno" ]] 
then 
	echo "Getting gene annotations..."
  resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
	minsample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)

  while read cmmd
  do
    echo $cmmd
    eval $cmmd

    sleep 1
  done < $PROJECT/Results/CMH/gene_anno_commands_resolution_"$resolution"_min_sample_"$minsample".txt
  exit
fi

##################################################################################
#####                          run model_case_sample.R   ######
##################################################################################
if [[ $1 = "caseinModel" ]]
then
  case_models=$PROJECT/Scripts/model_case_sample.R
  if [[ $2 = "--help" ]]
  then
    echo "Creating table with models..."
    # Learn about the input options
    $Rscript_var $case_models --help
  else
    qsub -wd $PROJECT -q centos78.q -o $PROJECT/Results/model_case_stdout.log -e $PROJECT/Results/model_case_stderr.log -pe smp 1 -V Scripts/run_collapsing_lclust.sh "caseinModel" $PROJECT
  fi  
  exit
fi

##################################################################################
# Run NHC. https://github.com/casanova-lab/NHC. Run on dev1:dev4
##################################################################################
if [[ $1 = "runNHC" ]]
then 
  model=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.model)
  cluster=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.cluster)
  add_param=$(cat ./Input/collapsing.yaml | shyaml get-value ADDITIONAL_PARAMETERS.runNHC)
  w=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.w)
  b=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.b)
  m=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.m)
  min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)

  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $run_NHC_var --help 
  else
    echo "Running NHC..."
    
    # run script
    $Rscript_var $run_NHC_var --model $model --resolution_var $resolution --cluster $cluster --mixed_ancestry_logical FALSE --m $m --b $b --w $w $add_param
  fi
  exit

fi

if [[ $1 = "runNHCMixedAncestry" ]]
then

  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var $run_NHC_var --help
  else
    model=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.model)
    w=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.w)
    b=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.b)
    m=$(cat ./Input/collapsing.yaml | shyaml get-value NHC_VARIABLE.m)
    mixed_ancestry_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.mixed_ancestry_cluster)
    max_ratio=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.max_ratio)
    min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
  
  
    # run script
    $Rscript_var $run_NHC_var --model $model --resolution_var $resolution --cluster $mixed_ancestry_cluster --min_sample $min_sample --max_ratio $max_ratio --mixed_ancestry_logical TRUE --b $b --w $w
  fi
  exit

fi
##################################################################################
# THE FOLLOWING CODE IS FOR TESTING NEW FEATURES. DO NOT USE UNLESS DEVELOPING...
##################################################################################
##################################################################################
# Find optimal LOEUF threshold
##################################################################################
if [[ $1 = "findThreshold" ]]
then
  echo "Finding optimal threshold..."
  if [[ $2 = "--help" ]]
  then
    # Learn about the input options
    $Rscript_var create_digenic_var=$PROJECT/Scripts/digenic_create_model.R
    $unbiased_threshold_var --help
  else
    cores=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.cores)
    qsub -wd $PROJECT -o $PROJECT/Results/stdout.log -e $PROJECT/Results/stderr.log -q centos78.q -pe smp $cores -V Scripts/run_collapsing_lclust.sh "findThreshold" $PROJECT
  fi

  exit
fi


##################################################################################
#####   The following steps are used for Regenie             ######
##################################################################################

##################################################################################
#####   Regenie VCF run on qs1   
#####   This should be run after Coverage2 ("clusterCoverage")
##################################################################################
## for these we may want to add --include-revel and --include-mtr to add as filters later on.
max_ratio=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.max_ratio)
mixed_ancestry_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.mixed_ancestry_cluster)
resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
regenie_folder="$PROJECT"/Results/Regenie/"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"

if [[ $1 = "RegListVCFMixedAncestry" ]]
then
  echo "In RegListVCFMixedAncestry"

  
  sleep 1
  samples=$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/*_existing.sample.txt
    
  chr="X"

  $atav --list-vcf --region $chr --gnomad-exome-filter-pass --gnomad-genome-filter-pass --gene-boundaries $PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/*site.clean.txt $QC --exclude-multiallelic-variant --sample $samples --out $regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/
  sleep 5

  for chr in {1..22}  # X done
    do

    $atav --list-vcf --region $chr --gnomad-exome-filter-pass --gnomad-genome-filter-pass --gene-boundaries $PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/*site.clean.txt $QC --exclude-multiallelic-variant --sample $samples --out $regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/
    sleep 5
  done
fi

# if [[ $1 = "RegListVCFMixedAncestry_VEP" ]]
# then
#   echo "In RegListVCFMixedAncestry_VEP"
# 
#   resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
# 
#   min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
#   
#   mixed_ancestry_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.mixed_ancestry_cluster)
#   
#   sleep 1
#   samples=$regenie_folder/samples.txt
#     
#   chr="21"
# 
#   # $atav --list-vcf --region $chr --gnomad-exome-filter-pass --gnomad-genome-filter-pass --gene-boundaries $PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/*site.clean.txt $QC --exclude-multiallelic-variant --sample $samples --out $regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/
#   
#   vcf=$regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/*.vcf.gz
#   echo $vcf
#   
# source /nfs/goldstein/software/centos7/vep_94.5_2020/vep_94.5-ENV.sh
# export PATH=/nfs/goldstein/software/vep_94.5:$PATH
# LD_LIBRARY_PATH=/nfs/goldstein/software/vep_94.5:$LD_LIBRARY_PATH
# source /usr/local/igm/non-atav-tools/htslib-1.10-x86_64/htslib-1.10-x86_64-ENV.sh
# source /nfs/goldstein/software/centos7/perl-5.32-x86_64/perl-5.32-ENV.sh
# 
#  vep -i $vcf --cache --offline --no_stats\
#     -dir /nfs/goldstein/software/vep_94.5/Cache \
#     --format vcf --vcf \
#     --assembly GRCh37 \
#     --fields "Allele,Consequence,Feature_type,Feature,SIFT,REVEL" --sift b \
#     -o $regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/Chr"$chr".annotated.vcf
# fi

if [[ $1 = "RegListVCFclustered" ]]
then
  mkdir -p $PROJECT/Regenie/
  touch $PROJECT/Regenie/samples.txt

  res_var=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution | sed 's/_/./')
  resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)

  min_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_cluster)
  max_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.max_cluster)
  sleep 1

  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do
    samples="$PROJECT/Results/KinshipFlashPCA/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
    cat $samples >> $PROJECT/Regenie/samples.txt
  done
  for (( i=$min_cluster; i<=$max_cluster; i++ ))
  do
    samples="$PROJECT/Results/KinshipFlashPCA/flashPCA_lclustering_res_"${resolution}"_cluster_"${i}"_sample.txt"
    
    chr="X"

    $atav --list-vcf --region $chr --gnomad-exome-filter-pass --gnomad-genome-filter-pass --gene-boundaries $PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}"_07/*site.clean.txt $QC --exclude-multiallelic-variant --sample $samples --out $PROJECT/Regenie/allQCBAFilter/allQCBAFilterChr${chr}_cluster_${i}/
    sleep 5

    for chr in {1..22}  # X done
    do

    $atav --list-vcf --region $chr --gnomad-exome-filter-pass --gnomad-genome-filter-pass --gene-boundaries $PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${i}"_07/*site.clean.txt $QC --exclude-multiallelic-variant --sample $samples --out $PROJECT/Regenie/allQCBAFilter/allQCBAFilterChr${chr}_cluster_${i}/
    sleep 5
    done
  done
  exit
fi

two_core_list="RegPLINKMergeMixedAncestry_1 RegPLINKMergeMixedAncestry_2 RegCovPheMixedAncestry_1 RegCovPheMixedAncestry_2 RegCovPheMixedAncestry_3 RegCovPheMixedAncestry_split RegStep1MixedAncestryBT RegStep2VCtests_QQManhattan RegListVCFMixedAncestry_create_AAF RegStep2VCtestsBT_ch RegStep2VCtestsBT_max RegStep2VCtestsBT_sum RegStep2ExWASBT"

if exists_in_list "$two_core_list" " " $1; then
  echo "$1 is a two cpu argument"
  echo "Running $1..."
  qsub -wd $PROJECT -o $PROJECT/Results/"$time_stamp"_stdout_$1.log -e $PROJECT/Results/"$time_stamp"_stderr_$1.log -q centos78.q -pe smp 2 -V Scripts/run_collapsing_lclust.sh $1 $PROJECT
  exit
else
  echo "more cpus"
fi


if [[ $1 = "RegStep2VCtestsBT" ]]
then
  $PROJECT/Scripts/collapsing_lclust.sh RegStep2VCtestsBT_ch
  sleep 2
  $PROJECT/Scripts/collapsing_lclust.sh RegStep2VCtestsBT_max
  sleep 2
  $PROJECT/Scripts/collapsing_lclust.sh RegStep2VCtestsBT_sum
fi
##################################################################################
#####   Regenie PLINK filtering/merge  run on qs1             ######
##################################################################################

##################################################################################
#####   Regenie Generate Covariate and phenotype files run on dev1-4        ######
##################################################################################

if [[ $1 = "RegListVCFMixedAncestry_genotype" ]]
then
  mkdir -p $regenie_folder
  mkdir $regenie_folder/allQCBAFilter
  mkdir $regenie_folder/allQCBAFilter/genotypes
  touch $regenie_folder/samples.txt

  samples="$PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/*existing.sample.txt"
  # cp $samples $regenie_folder/samples.txt
  # samples=$regenie_folder/samples.txt
  $atav --list-var-anno --out $regenie_folder/allQCBAFilter/genotypes --variant $regenie_folder/allQCBAFilter/allQCBAFilterVariantMap_just_id.txt --include-mtr --include-revel --include-rvis --include-sub-rvis --include-loftee --include-limbr --include-trap --include-primate-ai --include-pext --include-mpc --gzip
  exit
fi


# to adjust phenotypes, open up phenotypesBT.txt. All column, 1 case, 0 controls, NA will skip. option 1 is : "One Hot" encoding (column for each type of category, either 1 or 0 for each subtype) If doing smaller subset, do NAs. Add new column with refractory as 1, control as 0, non-refractory as NA. Alternatively, build a column with "types" where factors. (Refractory 1, non-refracotry 2, control 0). Regenie has feature where it will split it for you  

# if [[ $1 = "RegCovPheMixedAncestry_split" ]]


##################################################################################
#####   Regenie step 1 run on qs1 or dev1-4                                  #####
#####   Currently set up for binary traits, but can do quantitative traits as well #####
##################################################################################
# here, does whole genome regression. tries to account for ancestry, include covariates. output is not a statistical output. Just a setup for step 2. output files are not interpretable.
# if [[ $1 = "RegStep1MixedAncestryBT" ]]

if [[ $1 = "RegStep1MixedAncestryQT" ]]
then
  PROJECT="$PROJECT"
  mkdir $PROJECT/Results/Regenie/Step1New
  mkdir $PROJECT/Results/Regenie/Step1New/MinMAF1pct
  mkdir $PROJECT/Results/Regenie/Step1New/MinMAF1pct/QT

  step1=$PROJECT/Results/Regenie/Step1New/MinMAF1pct/QT
  cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
  pheno_file=$PROJECT/Results/Regenie/CovPheno/phenotypesQT.txt

  $regenie \
    --step 1 \
    --bed $PROJECT/Results/Regenie/allQCBAFilter/allQCBAFilter \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --qt --lowmem \
    --extract $PROJECT/Results/Regenie/allQCBAFilter/snps.pruned.QT.snplist \
    --bsize 1000 --loocv \
    --lowmem-prefix $step1/ \
    --out $step1/regenieRes
    # --exclude $PROJECT/Results/Regenie/allQCBAFilter/QT.exclude  \
  exit
fi

##################################################################################
#####   Regenie Step 2 approximate fast firth ExWAS run on dev1-4           ######
#####   This will run for binary traits outline in the pheno file          ######
##################################################################################



if [[ $1 = "RegStep2ExWASQT" ]]
then
bedfolder=$PROJECT/Results/Regenie/allQCBAFilter
outputFolder=$PROJECT/Results/Regenie/
bedfile=$bedfolder/allQCBAFilter.hwe.pruned

out=regenieRes.hwe.pruned
cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
pheno_file=$PROJECT/Results/Regenie/CovPheno/phenotypesQT.txt
step1=$PROJECT/Results/Regenie/Step1New/MinMAF1pct/QT
mkdir $outputFolder/ExWAS
mkdir $outputFolder/ExWAS/QT
outputFolder=$PROJECT/Results/Regenie/ExWAS/QT

minMac=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.minMac)
  echo "Minimum Minor Allele Count for regenie $minMac"

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --qt \
  --bsize 200 \
  --loocv \
  --pThresh 0.05 \
  --minMAC $minMac \
  --pred $step1/regenieRes_pred.list \
  --out $outputFolder/$out \
  --htp $outputFolder/$out
exit
fi

if [[ $1 = "RegStep2ExWASBT_SPA" ]]
then
bedfolder=$PROJECT/Results/Regenie/allQCBAFilter
outputFolder=$PROJECT/Results/Regenie/
bedfile=$bedfolder/allQCBAFilter.hwe.pruned

out=regenieRes.hwe.pruned_SPA
cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
pheno_file=$PROJECT/Results/Regenie/CovPheno/phenotypesBT.txt
step1=$PROJECT/Results/Regenie/Step1New/MinMAF1pct/BT
mkdir $outputFolder/ExWAS_SPA
mkdir $outputFolder/ExWAS_SPA/BT
outputFolder=$PROJECT/Results/Regenie/ExWAS_SPA/BT

minMac=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.minMac)
  echo "Minimum Minor Allele Count for regenie $minMac"

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --bt \
  --bsize 200 \
  --loocv \
  --spa \
  --pThresh 0.05 \
  --minMAC $minMac \
  --pred $step1/regenieRes_pred.list \
  --out $outputFolder/$out \
  --htp $outputFolder/$out
exit
fi


if [[ $1 = "RegStep2ExWASQT_SPA" ]]
then
bedfolder=$PROJECT/Results/Regenie/allQCBAFilter
outputFolder=$PROJECT/Results/Regenie/
bedfile=$bedfolder/allQCBAFilter.hwe.pruned

out=regenieRes.hwe.pruned_SPA
cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
pheno_file=$PROJECT/Results/Regenie/CovPheno/phenotypesQT.txt
step1=$PROJECT/Results/Regenie/Step1New/MinMAF1pct/QT
mkdir $outputFolder/ExWAS_SPA
mkdir $outputFolder/ExWAS_SPA/QT
outputFolder=$PROJECT/Results/Regenie/ExWAS_SPA/QT

minMac=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.minMac)
  echo "Minimum Minor Allele Count for regenie $minMac"

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --qt \
  --bsize 200 \
  --loocv \
  --spa \
  --pThresh 0.05 \
  --minMAC $minMac \
  --pred $step1/regenieRes_pred.list \
  --out $outputFolder/$out \
  --htp $outputFolder/$out
exit
fi


##################################################################################
#####   Regenie Step 2 approximate fast firth Burden run on dev1-4          ######
#####   This will run for binary traits outline in the pheno file           ######
####       DEPRECIATED by: "RegStep2VCtestsQT"                              ######
##################################################################################

# by gene association test.
if [[ $1 = "RegStep2BurdenQT" ]]
then
  PROJECT="$PROJECT"
  mkdir $PROJECT/Results/Regenie/Burden
  mkdir $PROJECT/Results/Regenie/Burden/QT
  outputFolder=$PROJECT/Results/Regenie/allQCBAFilter
  
  bedfolder=$PROJECT/Results/Regenie/allQCBAFilter
  bedfile=$bedfolder/allQCBAFilter.hwe.pruned
  
  out=regenieRes.hwe.pruned
  cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
  pheno_file=$PROJECT/Results/Regenie/CovPheno/phenotypesQTsc.txt
  step1=$PROJECT/Results/Regenie/Step1New/MinMAF1pct/QT
  outputFolder=$PROJECT/Results/Regenie/Burden/QT
  
  aaf="0.1,0.05,0.01"
  minMac=1
  
  mkdir $outputFolder/comphetAFcc${minMac}
  mkdir $outputFolder/maxAFcc${minMac}
  mkdir $outputFolder/sumAFcc${minMac}
  
  $regenie \
    --step 2 \
    --bed $bedfile \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --qt \
    --bsize 200 \
    --firth --approx \
    --pred $step1/regenieRes_pred.list \
    --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
    --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
    --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
    --aaf-bins $aaf \
    --minMAC $minMac \
    --build-mask "comphet" \
    --htp "htp" \
    --out $outputFolder/comphetAFcc${minMac}/$out
  
  $regenie \
    --step 2 \
    --bed $bedfile \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --qt \
    --bsize 200 \
    --firth --approx \
    --pred $step1/regenieRes_pred.list \
    --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
    --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
    --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
    --aaf-bins $aaf \
    --minMAC $minMac \
    --build-mask "max" \
    --htp "htp" \
    --out $outputFolder/maxAFcc${minMac}/$out
  
  $regenie \
    --step 2 \
    --bed $bedfile \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --qt \
    --bsize 200 \
    --firth --approx \
    --pred $step1/regenieRes_pred.list \
    --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
    --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
    --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
    --aaf-bins $aaf \
    --minMAC $minMac \
    --build-mask "sum" \
    --htp "htp" \
    --out $outputFolder/sumAFcc${minMac}/$out
  
  cd $PROJECT
  $Rscript_var $reg_results --aafs=$aaf
  
  echo "Regnie step 2 burden completed."
  exit
fi

##################################################################################
#####   Regenie Step 2 approximate fast firth Burden run on dev1-4          ######
#####   This will run for binary traits outline in the pheno file           ######
#####   IMPROTANT: This also runs a burden test as part of the VC testing   ######
#####        Therefore no need to run "RegStep2BurdenQT" on its own         ######
##################################################################################


# if [[ $1 = "RegStep2VCtestsBT" ]]

if [[ $1 = "RegStep2VCtestsQT" ]]
then
PROJECT="$PROJECT"
mkdir $PROJECT/Results/Regenie/VC
mkdir $PROJECT/Results/Regenie/VC/QT
outputFolder=$PROJECT/Results/Regenie/allQCBAFilter

bedfolder=$PROJECT/Results/Regenie/allQCBAFilter
bedfile=$bedfolder/allQCBAFilter.hwe.pruned

out=regenieRes.hwe.pruned
cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
pheno_file=$PROJECT/Results/Regenie/CovPheno/phenotypesQT.txt
step1=$PROJECT/Results/Regenie/Step1New/MinMAF1pct/QT
outputFolder=$PROJECT/Results/Regenie/VC/QT

aaf=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.aafs)
  echo "Max Minor Allele Fractions for regenie $aaf"
minMac=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.minMac)
  echo "Minimum Minor Allele Count for regenie $minMac"

mkdir $outputFolder/comphetAFcc${minMac}
mkdir $outputFolder/maxAFcc${minMac}
mkdir $outputFolder/sumAFcc${minMac}

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --qt \
  --bsize 200 \
  --firth --approx \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "comphet" \
  --htp "htp" \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/comphetAFcc${minMac}/$out

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --qt \
  --bsize 200 \
  --firth --approx \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "max" \
  --htp "htp" \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/maxAFcc${minMac}/$out

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --qt \
  --bsize 200 \
  --firth --approx \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "sum" \
  --htp "htp" \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/sumAFcc${minMac}/$out

echo "Regnie step 2 VC testing completed for QT."
exit
fi

if [[ $1 = "RegStep2VCtestsBT_SPA" ]]
then
PROJECT="$PROJECT"
mkdir $PROJECT/Results/Regenie/VC_SPA
mkdir $PROJECT/Results/Regenie/VC_SPA/BT
outputFolder=$PROJECT/Regenie/allQCBAFilter

bedfolder=$PROJECT/Results/Regenie/allQCBAFilter
bedfile=$bedfolder/allQCBAFilter.hwe.pruned

out=regenieRes.hwe.pruned_SPA
cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
pheno_file=$PROJECT/Results/Regenie/CovPheno/phenotypesBT.txt
step1=$PROJECT/Results/Regenie/Step1New/MinMAF1pct/BT
outputFolder=$PROJECT/Results/Regenie/VC_SPA/BT

aaf=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.aafs)
  echo "Max Minor Allele Fractions for regenie $aaf"
minMac=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.minMac)
  echo "Minimum Minor Allele Count for regenie $minMac"

mkdir $outputFolder/comphetAFcc${minMac}
mkdir $outputFolder/maxAFcc${minMac}
mkdir $outputFolder/sumAFcc${minMac}

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --bt \
  --bsize 200 \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "comphet" \
  --htp "htp" \
  --spa \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/comphetAFcc${minMac}/$out

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --bt \
  --bsize 200 \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "max" \
  --htp "htp" \
  --spa \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/maxAFcc${minMac}/$out

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --bt \
  --bsize 200 \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "sum" \
  --htp "htp" \
  --spa \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/sumAFcc${minMac}/$out

echo "Regnie step 2 VC testing completed for BT using SPA."
exit
fi

if [[ $1 = "RegStep2VCtestsQT_SPA" ]]
then
PROJECT="$PROJECT"
mkdir $PROJECT/Results/Regenie/VC_SPA
mkdir $PROJECT/Results/Regenie/VC_SPA/QT
outputFolder=$PROJECT/Results/Regenie/allQCBAFilter

bedfolder=$PROJECT/Results/Regenie/allQCBAFilter
bedfile=$bedfolder/allQCBAFilter.hwe.pruned

out=regenieRes.hwe.pruned_SPA
cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
pheno_file=$PROJECT/Results/Regenie/CovPheno/phenotypesQT.txt
step1=$PROJECT/Results/Regenie/Step1New/MinMAF1pct/QT
outputFolder=$PROJECT/Results/Regenie/VC_SPA/QT

aaf=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.aafs)
  echo "Max Minor Allele Fractions for regenie $aaf"
minMac=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.minMac)
  echo "Minimum Minor Allele Count for regenie $minMac"

mkdir $outputFolder/comphetAFcc${minMac}
mkdir $outputFolder/maxAFcc${minMac}
mkdir $outputFolder/sumAFcc${minMac}

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --qt \
  --bsize 200 \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "comphet" \
  --htp "htp" \
  --spa \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/comphetAFcc${minMac}/$out

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --qt \
  --bsize 200 \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "max" \
  --htp "htp" \
  --spa \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/maxAFcc${minMac}/$out

$regenie \
  --step 2 \
  --bed $bedfile \
  --covarFile $cov_file \
  --phenoFile $pheno_file \
  --qt \
  --bsize 200 \
  --pred $step1/regenieRes_pred.list \
  --anno-file $PROJECT/Results/Regenie/Anno/AnnoMerged.txt \
  --set-list $PROJECT/Results/Regenie/Anno/listFile.txt \
  --mask-def $PROJECT/Results/Regenie/Anno/mask.mask \
  --aaf-bins $aaf \
  --minMAC $minMac \
  --build-mask "sum" \
  --htp "htp" \
  --spa \
  --vc-tests skato-acat,acato-full \
  --out $outputFolder/sumAFcc${minMac}/$out

echo "Regnie step 2 VC testing completed for QT using SPA."
exit
fi

if [[ $1 = "RegStep2VCtests_QQManhattan_SPA" ]]
then
echo "Now splitting analysis by AF and Mask"

  aaf=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.aafs)
  echo "Max Minor Allele Fractions for regenie $aaf"
  proj_dir="$PWD"

 $Rscript_var $reg_results_SPA --aafs=$aaf # --project_dir=$proj_dir
exit
fi
