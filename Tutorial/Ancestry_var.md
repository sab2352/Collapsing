# Ancestry Var

#### The Ancestry_var filter is to select the minimum ancestry percentage for the Ancestry_filter option.

This filter is found within the <b> cohortSelection.R </b> script. 

#### The default is 0.75 if it is not set. This means that the samples need to have at least 75% match of ancestry. 

```
Usage: 

  cohortSelection.R --case_list_path=<case_list_path> [--ctrl_list_path=<ctrl_list_path>] [--sample_path_w_filename=<sample_path_w_filename>] [--case_group=<case_group>] [--exclude_list=<exclude_list>] [--BroadPhenotype_include=<BroadPhenotype_include>] [--acceptable_capture_kits=<acceptable_capture_kits>] [--DetailedPhenotype_excl_regex=<DetailedPhenotype_excl_regex>] [--CCDSBasesCov10X_var=<CCDSBasesCov10X_var>] [--DBSNPOverlapSNVs_var=<DBSNPOverlapSNVs_var>] [--DBSNPOverlapIndels_var=<DBSNPOverlapIndels_var>] [--ContaminationPercent_var=<ContaminationPercent_var>] [--Ancestry_var=<Ancestry_var>] [--Ancestry_filter=<Ancestry_filter>]

  
  Options:
  -h --help
  --case_list_path=<case_list_path> name of text file without header with list of cases to be included. Must be located in Input folder at top level of repo
  --ctrl_list_path=<ctrl_list_path> path to text file without header with list of ctrls to be included [default: NA]  
  --sample_path_w_filename=<sample_path_w_filename> path to csv file downloaded from https://sequence.igm.cumc.columbia.edu/ [default: /nfs/goldstein/software/atav_home/data/sample/igm_lims_sample_master_list.tsv]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --exclude_list=<exclude_list> specific samples to exclude [default: NA]  
  --BroadPhenotype_include=<BroadPhenotype_include> broad phenotypes to include in controls. Must uses _ instead of space. [default: healthy_family_member,control,control_mild_neuropsychiatric_disease]
  --acceptable_capture_kits=<acceptable_capture_kits> acceptable capture kits [default: NexteraRapidCapture,65MB,Roche,RocheV2,IDTERPv1,IDTERPv1Plus,IDTxGEN,IDTERPv1mtDNA,AgilentCRE,AgilentV4,Agilentv5,AgilentV5,AgilentV5UTR,AgilentV6,MedExome,N/A,Genome_v1,IDTERPv2]
  --DetailedPhenotype_excl_regex=<DetailedPhenotype_excl_regex> regexp string to exclude things in detailed phenotype [default: NA] 
  --CCDSBasesCov10X_var=<CCDSBasesCov10X_var> [default: 90] samples must have at least this % of the CCDS region covered by 10x
  --DBSNPOverlapSNVs_var=<DBSNPOverlapSNVs_var> [default: 0.85] must have at least this % overlap of SNVs with DBSNP
  --DBSNPOverlapIndels_var=<DBSNPOverlapIndels_var> [default: 0.80] must have at least this % overlap of indels with DBSNP
  --ContaminationPercent_var=<ContaminationPercent_var> [default: 2] upper limit to the contamination percentage for samples to include. Ignored if coverage >200x at 15% is used for these
  --debug=<debug> [default: FALSE]
  --Ancestry_filter=<Ancestry_filter> [default: NA] Predicted ancestry you wish to filter on
  --Ancestry_var=<Ancestry_var> [default: 0.75] Min ancestry percentage 

```
