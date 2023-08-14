## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>


# Creating Cohort

OPTION 1 (Default):
* The default sample_path_w_filename will be a continuously updated master sample file created by Nick. No action from user is needed. The only issue with this is that LOFTEE scores are only added monthly to LOFs new to atav
* https://redmine.igm.cumc.columbia.edu/issues/7624

OPTION 2:
* Go to https://sequence.igm.cumc.columbia.edu/search.php?action=searchSample
* Type "in dragen" into Current Prep Status Search and download results. 
* To download the result on a mac, hold option and click the download button. Add sample_path_w_filename argument in Yaml


**Runs on QS1**


Command: `$ ./Scripts/collapsing_lclust.sh cohortSelection`


```
Usage: 
  cohortSelection.R --case_list_path=<case_list_path> [--ctrl_list_path=<ctrl_list_path>] [--sample_path_w_filename=<sample_path_w_filename>] [--case_group=<case_group>] [--exclude_list=<exclude_list>] [--BroadPhenotype_include=<BroadPhenotype_include>] [--acceptable_exome_kits=<acceptable_exome_kits>] [--DetailedPhenotype_excl_regex=<DetailedPhenotype_excl_regex>] [--CCDSBasesCov10X_var=<CCDSBasesCov10X_var>] [--DBSNPOverlapSNVs_var=<DBSNPOverlapSNVs_var>] [--DBSNPOverlapIndels_var=<DBSNPOverlapIndels_var>]
  
  Options:
  -h --help
  --case_list_path=<case_list_path> path to text file without header with list of cases to be included
  --ctrl_list_path=<ctrl_list_path> path to text file without header with list of ctrls to be included [default: NA]  
  --sample_path_w_filename=<sample_path_w_filename> path to csv file downloaded from https://sequence.igm.cumc.columbia.edu/ [default: /nfs/goldstein/software/atav_home/data/sample/igm_lims_sample_master_list.tsv]
  --case_group=<case_group> label for analysis. example might "picu" or "IPF" [default: temp_case]
  --exclude_list=<exclude_list> specific samples to exclude [default: NA]  
  --BroadPhenotype_include=<BroadPhenotype_include> broad phenotypes to include in controls. Must uses _ instead of space. [default: healthy_family_member,control,control_mild_neuropsychiatric_disease]
  --acceptable_exome_kits=<acceptable_exome_kits> acceptable exome capture kits [default: NexteraRapidCapture,65MB,Roche,RocheV2,IDTERPv1,IDTERPv1Plus,IDTxGEN,IDTERPv1mtDNA,AgilentCRE,AgilentV4,Agilentv5,AgilentV5,AgilentV5UTR,AgilentV6,MedExome,N/A,Genome_v1]
  --DetailedPhenotype_excl_regex=<DetailedPhenotype_excl_regex> regexp string to exclude things in detailed phenotype [default: abcdefg] 
  --CCDSBasesCov10X_var=<CCDSBasesCov10X_var> [default: 90] samples must have at least this % of the CCDS region covered by 10x
  --DBSNPOverlapSNVs_var=<DBSNPOverlapSNVs_var> [default: 0.85] must have at least this % overlap of SNVs with DBSNP
  --DBSNPOverlapIndels_var=<DBSNPOverlapIndels_var> [default: 0.80] must have at least this % overlap of indels with DBSNP  
```
### Yaml Variables Used:
```
USER_VARIABLE:
  case_path: /PATH/case_list.txt
  case_group: project_name
ADDITIONAL_PARAMETERS:
  cohortSelection: ""
```
* case_path: provides path to case list
* case_group: project name, added to output file(s) name
* cohortSelection: enter any additional parameters you wish to specify between the quotes
  * should be look like cohortSelection: "--DBSNPOverlapSNVs_var 0.85 --BroadPhenotype_include healthy_family_member,control"

Input
- Will use the case list txt file

Output
- /Data
    - .csv files that contain case list and controls
    - .ped.txt file that contain case and control sample names
    - .txt file that lists cases missing from dragon
