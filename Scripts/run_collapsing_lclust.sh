#$ -S /bin/sh   
#$ -M myuni@cumc.columbia.edu
#$ -m bea
#$ -N syn

# shyaml requires this version of python. This can be sourced on qs1
source /nfs/goldstein/software/centos7/python-3.9.7-x86_64_shared/python3.9.7-ENV.sh
# this makes sure that we're all using the same version of Ron the server.
source /nfs/goldstein/software/centos7/R-4.1.0_with_gcc_10-x86_64/R-4.1.0_with_gcc_10-x86_64-ENV.sh
plink2="/usr/local/igm/non-atav-tools/plink2_20220814/plink2"
regenie="/nfs/goldstein/software/centos7/regenie-3.1.3/regenie_v3.1.3.gz_x86_64_Centos7_mkl"

echo "in run_collapsing_lclust.sh"

PROJECT=$2
case_group=$(cat ./Input/collapsing.yaml | shyaml get-value USER_VARIABLE.case_group)
echo "case_group is $case_group"

if [[ $1 = "cohortSelection" ]]
then

  case_path=$(cat ./Input/collapsing.yaml | shyaml get-value USER_VARIABLE.case_path)
  echo "case_path path is "
  echo $case_path
  echo "Creating cohort ..."
  add_param=$(cat ./Input/collapsing.yaml | shyaml get-value ADDITIONAL_PARAMETERS.cohortSelection)

  # run cohort selection
  Rscript $PROJECT/Scripts/cohortSelection.R --case_list_path $case_path --case_group $case_group $add_param
  exit
fi

if [[ $1 = "clustering" ]]
then
  resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution  | sed 's/_/./')
  echo "Creating clusters..."
  echo "Resolution"
  echo $resolution
  # run cluster creation
  Rscript $PROJECT/Scripts/lclust_Flash.R --resolution_var $resolution --case_group $case_group
  exit
fi

if [[ $1 = "cmh" ]]
then
	resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
	echo "Resolution is $resolution"
	minsample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
	echo "Min Sample is $minsample"
	models_to_exclude=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.models_to_exclude)
	echo "models_to_exclude is $models_to_exclude"

	# run script
	Rscript $PROJECT/Scripts/cmh_lclust.R --resolution_var $resolution --min_sample $minsample --models_to_exclude $models_to_exclude
  exit
fi


if [[ $1 = "cmhPerm" ]]
then

	cmh_perm=$PROJECT/Scripts/cmh_permute_lclust_short_script.R
	resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
	minsample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
	model=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.perm_models)
	echo "Models for cmh_permute are $model"
	permstart=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.permstart)
	permend=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.permend)
	cores=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.cores)

	Rscript $cmh_perm --project $PROJECT/Results  --res $resolution --minsample $minsample --model $model --permstart $permstart --permend $permend --cores $cores
  exit
fi

if [[ $1 = "create_digenic" ]]
then
	create_digenic_var=$PROJECT/Scripts/digenic_create_model.R
	resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)

	Rscript $create_digenic_var --resolution_var $resolution --model $3
    exit
fi

if [[ $1 = "findThreshold" ]]
then
	unbiased_threshold_var=$PROJECT/Scripts/unbiased_threshold.R
	resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
	model=$(cat ./Input/collapsing.yaml | shyaml get-value THRESHOLD_DETECTION.model)
	echo $model
	permend=$(cat ./Input/collapsing.yaml | shyaml get-value THRESHOLD_DETECTION.permend)
	echo $permend
	threshold_var=$(cat ./Input/collapsing.yaml | shyaml get-value THRESHOLD_DETECTION.threshold_var)
	echo $threshold_var
  	cores=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.cores)
	minsample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
  
  Rscript $unbiased_threshold_var --model $model --resolution_var $resolution --min_sample $minsample --threshold_var $threshold_var --nperm $permend --cores $cores --case_group $case_group

  exit
fi

if [[ $1 = "create_gene_set_model" ]]
then
  	create_gene_set=$PROJECT/Scripts/fp_create_gene_set.R
  	echo "Creating gene sets ..."
  	resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
  	gene_set_path=$(cat ./Input/collapsing.yaml | shyaml get-value GENE_SET.path)
    model=$(cat ./Input/collapsing.yaml | shyaml get-value GENE_SET.model)
    gene_set_name=$(cat ./Input/collapsing.yaml | shyaml get-value GENE_SET.name)

  	Rscript $create_gene_set --model $model --resolution $resolution --cluster $3 --gene_set_path $gene_set_path --gene_set_name $gene_set_name --ccds_genes_path $4


	exit
fi

if [[ $1 = "filter_loeuf" ]]
then
  	echo "LOEUF thresholding ..."
    resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
    model=$(cat ./Input/collapsing.yaml | shyaml get-value LOEUF_FILTER.model)
    loeuf=$(cat ./Input/collapsing.yaml | shyaml get-value LOEUF_FILTER.loeuf)

  	Rscript $PROJECT/Scripts/fp_filter_loeuf.R --model $model --resolution $resolution --cluster $3 --loeuf $loeuf


	exit
fi

if [[ $1 = "forestPlot_display" ]]
then
	echo "Creating forest plot..."
	add_param=$(cat ./Input/collapsing.yaml | shyaml get-value ADDITIONAL_PARAMETERS.forestPlot)
	path_to_fp_models=$(cat ./Input/collapsing.yaml | shyaml get-value FOREST_PLOT.path_to_fp_models)
	cores=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.cores)

	mkdir $PROJECT/Results/fp/
	Rscript $PROJECT/Scripts/fp_forest_plot.R --output_directory $PROJECT/Results/fp/ --path_to_fp_models $path_to_fp_models --cores 2 $add_param
  exit
fi

if [[ $1 = "TableWModels" ]]
then

    echo "Creating table with models..."
    add_param=$(cat ./Input/collapsing.yaml | shyaml get-value ADDITIONAL_PARAMETERS.TableWModels)
    model_name=$(cat ./Input/collapsing.yaml | shyaml get-value TABLE_W_MODELS.model_name)
    inheritance=$(cat ./Input/collapsing.yaml | shyaml get-value TABLE_W_MODELS.inheritance)
    effects=$(cat ./Input/collapsing.yaml | shyaml get-value TABLE_W_MODELS.effects)
    
    echo "$add_param $model_name $inheritance $effects"
    # run script
    Rscript $PROJECT/Scripts/TableWModels.R --model_name $model_name --inheritance $inheritance --effects $effects  $add_param --case_group $case_group
  exit
fi

if [[ $1 = "sampleTracking" ]]
then
    echo "Creating sample tracking log..."
    resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
	min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
    
    # run script
    Rscript $PROJECT/Scripts/table_sample_counts.R --resolution_var $resolution --min_sample $min_sample --case_group $case_group
  exit
fi

if [[ $1 = "qq" ]]
then
	cmh_qq_var=$PROJECT/Scripts/cmh_qq.R
	min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
	resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
	permend=$(cat ./Input/collapsing.yaml | shyaml get-value PERMUTATION_VARIABLE.permend)
	add_param=$(cat ./Input/collapsing.yaml | shyaml get-value ADDITIONAL_PARAMETERS.cmh_qq)
    echo "Creating qq plots ..."
    
    # run script
    Rscript $PROJECT/Scripts/cmh_qq.R --resolution_var $resolution --min_sample $min_sample --nperm_end $permend --case_group $case_group $add_param
  exit
fi

##################################################################################
#####                          combine cluster genotype files             ######
##################################################################################
if [[ $1 = "combine_cluster_genotype" ]]
then
  resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
  min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
  combine_ctrl_output=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.combine_ctrl_output)

  Rscript $PROJECT/Scripts/review_genotypes_from_model.R --resolution $resolution --min_sample $min_sample --case_group $case_group --combine_ctrl_output $combine_ctrl_output
  exit
fi

if [[ $1 = "summarizeRemovedSamples" ]]
then
	case=$(cat ./Input/collapsing.yaml | shyaml get-value USER_VARIABLE.case_path)
  
	python $PROJECT/Scripts/identify_removed_samples.py --caselist $PROJECT/Input/$case --ped_txt $PROJECT/Data/*.ped.txt --sample_file_check $PROJECT/Results/SampleFileCheck/*_existing.sample.txt --out $PROJECT/Results/SampleFileCheck/Check
  exit  
fi

if [[ $1 = "caseinModel" ]]
then
	min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
	resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
  echo "Identifying case samples present in models..."
    
  # run script
	Rscript $PROJECT/Scripts/model_case_sample.R --resolution_var $resolution --min_sample $min_sample 
  exit
fi
min_sample=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.min_sample)
max_ratio=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.max_ratio)
resolution=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.resolution)
if [[ $1 = "clusteringMixedAncestry" ]]
then
  echo "Create mixed ancestry cluster..."

  # run cluster creation
  Rscript $PROJECT/Scripts/lclust_Flash_combo_cluster.R --resolution_var $resolution --case_group $case_group --min_sample $min_sample --max_ratio $max_ratio
  exit
fi

########
# REGENIE
########
mixed_ancestry_cluster=$(cat ./Input/collapsing.yaml | shyaml get-value COLLAPSING_VARIABLE.mixed_ancestry_cluster)
regenie_folder="$PROJECT"/Results/Regenie/"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"
plink_mac=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.plink_mac)
pcs=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.pcs)  


# For this step... make bed files from VCF files... Merge BED files... extract common variants (> MAF 0.001 or 0.01 depending on code)... then removes variants that are out of HWE
plink="/nfs/goldstein/software/PLINK/PLINK_1.90_3.38/plink"
if [[ $1 = "RegPLINKMergeMixedAncestry_1" ]]
then
  # Create bed file from vcf for chromosome x. See https://zzz.bwh.harvard.edu/plink/binary.shtml for description of .bed file format
  chr="X"
  
  file=$regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/*.vcf.gz
  
  file_new=$regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/allQCBAFilterChr$chr
  
  $plink --vcf $file --make-bed --out $file_new

  # Create bed file from vcf for autosomes. See https://zzz.bwh.harvard.edu/plink/binary.shtml for description of .bed file format
  for chr in {1..22} 
  do
  
  file=$regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/*.vcf.gz
  
  file_new=$regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/allQCBAFilterChr$chr
  
  $plink --vcf $file --make-bed --out $file_new
  
  done
  
  # Find the .bim files. bim files are extended map files, which also includes the names of the alleles: (chromosome, SNP, cM, base-position, allele 1, allele 2). see https://zzz.bwh.harvard.edu/plink/binary.shtml
  cd $regenie_folder/allQCBAFilter
  
  find "$(pwd -P)" -type f -name "*QCBAFilterChr*.bim" > $regenie_folder/allQCBAFilter/ForMerge.list
  
  
  sed -i 's/.bim$//g' $regenie_folder/allQCBAFilter/ForMerge.list
  #GOTO_1
  # merge all files in one. ForMerge.list is a list of folders, not files.
  # --merge-list <filename> This allows you to merge more than one fileset to the reference fileset. (Also, this can be used without a reference; in that case, the newly created fileset is then treated as the reference by most other PLINK operations.) The parameter must be the name of a text file specifying one fileset per line.
  outputFolder=$regenie_folder/allQCBAFilter
  $plink --merge-list $regenie_folder/allQCBAFilter/ForMerge.list --out $outputFolder/allQCBAFilter
  
  chr=X
# This command searches for lines in a file that do not contain "##", selects columns 3 and 8 from those lines, and saves the extracted data into a file named allQCBAFilterChr${chr}variants.txt. In this case, row 3 is the variant ID and row 8 is the info. You end up with a text file with just those two columns. 
  file=$regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/*allQCBAFilterChr$chr*.vcf.gz
  zgrep -v "##" $file|cut -f 3,8 > $regenie_folder/allQCBAFilter/allQCBAFilterChr${chr}variants.txt
  
  for chr in {1..22}
  do
    echo "exporting variants from chr ${chr}"
    file=$regenie_folder/allQCBAFilter/allQCBAFilterChr$chr/*allQCBAFilterChr$chr*.vcf.gz
    zgrep -v "##" $file|cut -f 3,8 > $regenie_folder/allQCBAFilter/allQCBAFilterChr${chr}variants.txt
  
  done
  echo "Done with $1"
  exit
fi

if [[ $1 = "RegListVCFMixedAncestry_create_AAF" ]]
then
  geno_file="$regenie_folder/allQCBAFilter/genotypes/*genotypes_annotations.csv"

  Rscript $PROJECT/Scripts/Regenie_make_aaf_file.R --case_group $case_group --regenie_folder $regenie_folder
  exit
fi

bfile=$regenie_folder/allQCBAFilter/allQCBAFilter
if [[ $1 = "RegPLINKMergeMixedAncestry_2" ]]
then
  # The --bfile flag causes the binary fileset plink.bed + plink.bim + plink.fam to be referenced. 
  # To exclude markers that failure the Hardy-Weinberg test at a specified significance threshold, use the option:
  # plink --file mydata --hwe 0.001
  # --hardy writes a list of genotype counts and Hardy-Weinberg equilibrium exact test statistics to plink.hwe. With the 'midp' modifier, a mid-p adjustment is applied 
  # --hwe's 'midp' modifier applies the mid-p adjustment described in Graffelman J, Moreno V (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium. The mid-p adjustment tends to bring the null rejection rate in line with the nominal p-value, and also reduces the filter's tendency to favor retention of variants with missing data. We recommend its use.
  # need discussion regarding HWE effect
  # --make-bed creates a new PLINK 1 binary fileset, after applying sample/variant filters and other operations below. For example, Generate binary_fileset.bed + .bim + .fam. Any samples/variants removed from the current analysis are also not present in this fileset. (This is the --make-bed step.)
  $plink --bfile $bfile \
    --hwe 1E-30 midp \
    --make-bed \
    --out $bfile.hwe.pruned 
  exit
fi

if [[ $1 = "RegCovPheMixedAncestry_1" ]]
then
  echo "Creating Regenie Covariate and Phenotype files..."

  # Learn about the input options
  # $Rscript_var $reg_covphe --help

  ## run script
  echo "Running regenieCovPheno.R"
  Rscript $PROJECT/Scripts/regenieCovPheno.R \
    --pcs $pcs \
    --pc_path "$PROJECT/Results/KinshipFlashPCA/PCAFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/" \
    --coverage_path $PROJECT/Results/Coverage/CoverageFlashLClust_res_"${resolution}"_cluster_"${mixed_ancestry_cluster}"_min_sample_"${min_sample}"_max_ratio_"${max_ratio}"_07_mixedCluster/ \
    --regenie_folder $regenie_folder
  echo "Finished regenieCovPheno.R"
  exit
fi

if [[ $1 = "RegCovPheMixedAncestry_2" ]]
then
  echo "Running AnnoMaker.sh"
  $PROJECT/Scripts/AnnoMaker.sh $regenie_folder
  echo "Finished AnnoMaker.sh"
  exit
fi
if [[ $1 = "RegCovPheMixedAncestry_3" ]]
then
  echo "Running regenieCovPheno_update.R"
  Rscript $PROJECT/Scripts/regenieCovPheno_update.R --regenie_folder $regenie_folder --case_group $case_group
  echo "Finished regenieCovPheno_update.R"
  exit
fi

if [[ $1 = "RegCovPheMixedAncestry_split" ]]
then
  echo "RegCovPheMixedAncestry_split..."
  echo "Running regenieCovPhenoSplit.R..."
  Rscript $PROJECT/Scripts/regenieCovPhenoSplit.R \
    --case_group $case_group \
    --pcs $pcs \
    --regenie_folder $regenie_folder
  echo "Finished regenieCovPhenoSplit.R..."
  
  ## Binary traits
  echo "Running plink2 binary trait"
  # --maf filters out all variants with allele frequency below the provided threshold (default 0.01), while --max-maf imposes an upper bound. Similarly, --mac and --max-mac impose lower and upper allele count bounds, respectively.
  # By default, unless the input is loaded with --no-sex1, samples with ambiguous sex have their phenotypes set to missing when analysis commands are run. Use --allow-no-sex to prevent this. (This setting is no longer ignored when --make-bed or --recode is present.)
  # This code filters out singletons. Josh is wondering whether we should be using --maf with either 0.01 or 0.001 in addition since the snp list passing that threshold which is generated preivously is never used.
  $plink2 --bfile $regenie_folder/allQCBAFilter/allQCBAFilter.hwe.pruned --allow-no-sex --pheno $regenie_folder/CovPheno_pc"$pcs"/phenotypesBT.txt --prune --mac $plink_mac --write-snplist --out $regenie_folder/allQCBAFilter/snps_pc"$pcs"_mac_"$plink_mac".pruned.BT
  ## Quantitative traits
  echo "Running plink2 quant trait"
  # $plink2 --bfile $regenie_folder/allQCBAFilter/allQCBAFilter.hwe.pruned --allow-no-sex --pheno $regenie_folder/CovPheno_pc"$pcs"/phenotypesQT.txt --prune --mac $plink_mac --write-snplist --out $regenie_folder/allQCBAFilter/snps_pc"$pcs"_mac_"$plink_mac".pruned.QT
  # naama version
  # $plink2 --bfile $PROJECT/Results/Regenie/allQCBAFilter/allQCBAFilter.hwe.pruned --allow-no-sex --1 --pheno $PROJECT/Results/Regenie/CovPheno/phenotypesBT.txt --maf 0.01 --mac 3 --write-snplist --out $PROJECT/Results/Regenie/allQCBAFilter/snps.pruned.BT
  exit
fi

step1=$regenie_folder/Step1New_pc"$pcs"_mac_"$plink_mac"/BT
mkdir -p $step1
if [[ $1 = "RegStep1MixedAncestryBT" ]]
then
  echo "RegStep1MixedAncestryBT..."
  
  # Running --step 1 --out foo will produce
    # A set of files containing genomic predictions for each phenotype from Step 1 (see Output section below).
    # A file called foo_pred.list listing the locations of the prediction files.


  # mkdir $regenie_folder/Step1New
  # mkdir $regenie_folder/Step1New/MinMAF1pct
  # mkdir $regenie_folder/Step1New/MinMAF1pct/BT
  # 
  # cov_file=$PROJECT/Results/Regenie/CovPheno/covariates.txt
  cov_file=$regenie_folder/CovPheno_pc"$pcs"/covariates.txt
  pheno_file=$regenie_folder/CovPheno_pc"$pcs"/phenotypesBT.txt

  # --lowmem flag to reduce memory usage by writing level 0 predictions to disk (details below). This is very useful if the number of traits is large (e.g. greater than 10)
  # --bt specify that traits are binary with 0=control,1=case,NA=missing
  # --extract	FILE	Optional	Inclusion file that lists IDs of variants to keep. I think these should be the common variants in HWE but I'm not sure if we're filtering them adequately for common variants.
  # --bsize size of the genotype blocks
  # --loocv	FLAG	Optional	flag to use leave-one out cross validation
  $regenie \
    --step 1 \
    --bed $regenie_folder/allQCBAFilter/allQCBAFilter \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --bt --lowmem \
    --extract $regenie_folder/allQCBAFilter/snps_pc"$pcs"_mac_"$plink_mac".pruned.BT.snplist \
    --bsize 1000 --loocv \
    --lowmem-prefix $step1/ \
    --out $step1/regenieRes
    # --exclude $PROJECT/Results/Regenie/allQCBAFilter/BT.exclude \
  exit
fi

vc_folder=$regenie_folder/VC_pc"$pcs"_mac_"$plink_mac"
aaf_file="$regenie_folder/allQCBAFilter/genotypes/aaf_file.txt"
bedfolder=$regenie_folder/allQCBAFilter
bedfile=$bedfolder/allQCBAFilter.hwe.pruned
cov_file=$regenie_folder/CovPheno_pc"$pcs"/covariates.txt
pheno_file=$regenie_folder/CovPheno_pc"$pcs"/phenotypesBT.txt
outputFolder=$vc_folder/BT
mkdir -p $outputFolder
out=regenieRes.hwe.pruned

aaf=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.aafs)
  echo "Max Minor Allele Fractions for regenie $aaf"
minMac=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.minMac)
  echo "Minimum Minor Allele Count for regenie $minMac"

if [[ $1 = "RegStep2VCtestsBT_ch" ]]
then
  echo "$1 ..."

  
  mkdir $outputFolder/comphetAFcc${minMac}
# --firth	FLAG	Optional	specify to use Firth likelihood ratio test (LRT) as fallback for p-values less than threshold
# --approx	FLAG	Optional	flag to use approximate Firth LRT for computational speedup (only works when option --firth is used)
# --pred	FILE	Optional	File containing predictions from Step 1 (see Overview). This is required for --step 2
# --anno-file	FILE	Required	File with variant annotations for each set
# --set-list	FILE	Required	File listing variant sets
# --mask-def	FILE	Required	File with mask definitions using the annotations defined in --anno-file
# --minMAC	FLOAT	Optional	flag to specify the minimum minor allele count (MAC) when testing variants [default is 5]. Variants with lower MAC are ignored. https://github.com/rgcgithub/regenie/issues/89 i think this is mask based
# Option --aaf-bins specifies the AAF upper bounds used to generate burden masks (AAF and not MAF [minor allele frequency] is used when deciding which variants go into a mask). By default, a mask based on singleton sites are always included. For example, --aaf-bins 0.01,0.05 will generate 3 burden masks for AAFs in [0,0.01], [0,0.05] and singletons.
# --build-mask	STRING	Optional	build masks using the maximum number of ALT alleles across sites ('max'; the default), or the sum of ALT alleles ('sum'), or thresholding the sum to 2 ('comphet')


  $regenie \
    --step 2 \
    --bed $bedfile \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --bt \
    --bsize 200 \
    --firth --approx \
    --pred $step1/regenieRes_pred.list \
    --anno-file $regenie_folder/Anno/AnnoMerged.txt \
    --set-list $regenie_folder/Anno/listFile.txt \
    --mask-def $regenie_folder/Anno/mask.mask \
    --aaf-bins $aaf \
    --aaf-file $aaf_file \
    --minMAC $minMac \
    --build-mask "comphet" \
    --write-mask-snplist \
    --htp "htp" \
    --vc-tests skato-acat,acato-full \
    --out $outputFolder/comphetAFcc${minMac}/$out
    
  echo "Completed $1"
  exit
fi

if [[ $1 = "RegStep2VCtestsBT_max" ]]
then
  echo "$1 ..."

  mkdir $outputFolder/maxAFcc${minMac}

  $regenie \
    --step 2 \
    --bed $bedfile \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --bt \
    --bsize 200 \
    --firth --approx \
    --pred $step1/regenieRes_pred.list \
    --anno-file $regenie_folder/Anno/AnnoMerged.txt \
    --set-list $regenie_folder/Anno/listFile.txt \
    --mask-def $regenie_folder/Anno/mask.mask \
    --aaf-bins $aaf \
    --aaf-file $aaf_file \
    --write-mask-snplist \
    --minMAC $minMac \
    --build-mask "max" \
    --htp "htp" \
    --vc-tests skato-acat,acato-full \
    --out $outputFolder/maxAFcc${minMac}/$out

  echo "Completed $1"
  exit
fi

if [[ $1 = "RegStep2VCtestsBT_sum" ]]
then
  echo "$1 ..."

  mkdir $outputFolder/sumAFcc${minMac}
  
  $regenie \
    --step 2 \
    --bed $bedfile \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --bt \
    --bsize 200 \
    --firth --approx \
    --pred $step1/regenieRes_pred.list \
    --anno-file $regenie_folder/Anno/AnnoMerged.txt \
    --set-list $regenie_folder/Anno/listFile.txt \
    --mask-def $regenie_folder/Anno/mask.mask \
    --aaf-file $aaf_file \
    --aaf-bins $aaf \
    --minMAC $minMac \
    --build-mask "sum" \
    --htp "htp" \
    --vc-tests skato-acat,acato-full \
    --out $outputFolder/sumAFcc${minMac}/$out
  
  echo "Completed $1"
  exit
fi

if [[ $1 = "RegStep2VCtests_QQManhattan" ]]
then
echo "Now splitting analysis by AF and Mask"

  aaf=$(cat ./Input/collapsing.yaml | shyaml get-value REGENIE_VARIABLE.aafs)
  echo "Max Minor Allele Fractions for regenie $aaf"
  
  Rscript $PROJECT/Scripts/regenieResults.R \
    --regenie_folder $regenie_folder \
    --vc_folder $vc_folder \
    --pcs $pcs \
    --case_group $case_group \
    --aafs=$aaf 
  exit
fi

# This will run single variant analysis
if [[ $1 = "RegStep2ExWASBT" ]]
then
  mkdir $vc_folder/ExWAS
  mkdir $vc_folder/ExWAS/BT

  $regenie \
    --step 2 \
    --bed $bedfile \
    --covarFile $cov_file \
    --phenoFile $pheno_file \
    --bt \
    --bsize 200 \
    --loocv \
    --aaf-file $aaf_file \
    --firth --approx \
    --pThresh 0.05 \
    --minMAC $minMac \
    --pred $step1/regenieRes_pred.list \
    --out $vc_folder/ExWAS/BT/$out \
    --htp $vc_folder/ExWAS/BT/$out
  exit
fi
