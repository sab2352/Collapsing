## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Run Regenie

We have included the scripts required to run regenie in this repository, using the ATAV outputs. \
For more information about regenie you can visit their [documentation](https://rgcgithub.github.io/regenie/) or [github](https://github.com/rgcgithub/regenie) pages.

The first step is to correctly set up the `Input/collapsing.yaml` file. \
Right now there are 2 regenie related variables in the yaml: The maximum alternate allele frequence cut off to use (aafs) and the minimum minor allele count for a gene (minMac). \
These are set at the defaults but can be changed. 

The first steps follow the standard collapsing steps including:
1. [Repository Setup](/Tutorial/Repository-Setup.md)
2. [Creating Cohort](/Tutorial/Creating-Cohort.md)
3. [First Coverage of CaseControl Cohort](/Tutorial/Coverage-of-CaseControl-Cohort.md)
4. [Principle Componenet Analysis and Removal of Related Samples](/Tutorial/Creation-of-Principle-Components-and-Removal-of-Related-Samples.md)

    At this point the script is set up to use a single multi-ancestry cluster. To do this I recommend doing regular clustering first. 
5. [Create Clusters](/Tutorial/Create-Clusters.md) 
Then you can set the `Input/collapsing.yaml` value COLLAPSING_VARIABLE:mixed_ancestry_cluster: to the number of clusters +1 (eg. if you create clusters 0 through 9, then `mixed_ancestry_cluster: 10`). Here you can also update the max_cluster value to the number of initial clusters created (eg. as per above this would be set as 9). We can also adjust the max_ratio and min_sample. Given regenie deals well with case:control imbalance, it is ok to set `max_ratio: NA`. 
6. [Create Mixed Ancestry Cluster](/Tutorial/Create-Mixed-Ancestry-Cluster.md). 
7. Run a PCA on the mixed-ancestry cluster using `./Scripts/collapsing_lclust.sh pcaMixedAncestry` run on qs1. These are used in the creation of the covariate file to account for ancestry differences.
8. Coverage in the Mixed Ancestry Cluster. This is done by calling `./Scripts/collapsing_lclust.sh clusterCoverageMixedAncestry` run on qs1

From Here we move into Regenie specific commands: 
1. Make VCF files from our mixed ancestry cluster. 
    - This is done with `./Scripts/collapsing_lclust.sh RegListVCFMixedAncestry ` run on the qs1 server.
    - The output of this is a single VCF for each chromosome, and it copies your sample list into the file: `Regenie/samples.txt`
    - Also runs a `--list-var-geno` which is used later for variant intolerance scores
2. Merge the VCF files.
    - This is done with `./Scripts/collapsing_lclust.sh RegPLINKMergeMixedAncestry ` run on the qs1 server.
    - If you want to change the Hardy-Weinburg filtering used, this is the step to change it via altering the code: `--hwe 1E-30 midp`.
    - Can adjust down if you want to `--hwe 1E-7 midp` which will filter more variants but has a lower risk of false positives from ancestry driven signals.
    - The output of this is a plink bed file and variant list file for each chromosome (in their respective chromosome folders), and in `Regenie/allQCBAFilter/allQCBAFilter` a merged plink bed file, a HWE pruined bed file, and a maf 0.001 filtered snplist file to use in regenie step 1. 
    - I found using a maf of 0.01 resulted in too few variants for step 1 when using exomes and hundreds of cases. If you have WGS or arrays + exomes or larger samples a maf of 0.01 may be more appropriate.
3. Make the covariate, phenotype, and variant type filter files:
    - This is done with `./Scripts/collapsing_lclust.sh RegCovPheMixedAncestry `
    - I would recommend running this with a nohup as it takes 3-4 hours to run and if the terminal loses connection **it will SILENTLY FAIL**. eg.  `nohup ./Scripts/collapsing_lclust.sh RegCovPheMixedAncestry &` (This advice is no longer necessary as now submitted to qs1)
    - You need to look at the elbow plot created from your PCA to determine the optimal PCs to include. Default covariates are sex and first 6 PCs. This can be increased or descreased by adjusting the `REGENIE_VARIABLE.pcs` value in the collapsing.yaml file to the number of PC dimensions to include.
    - The phenotype is currently set to run for binary traits defined by the V6 (6th) column of your fam file to compare cases with controls.
    - **You LIKELY will need to do some hands on coding of the phenotype file. It can accept any binary traist (coded 0=control,1=case,NA=excluded), quantitative traits (any numeric value=included, NA=excluded). Many phenotypes (100-1000s) can be analyzed in parallel.**
    - Regenie is however limited to looking at only binary or quantitative traits in a single analysis, so if we have both, we need to split them going forward.
    - This also created the masks and filters required for the burden and VC testing used later. This can be changed in the future
4. Splitting of Binary and Quantitative traits
    - This is done with `./Scripts/collapsing_lclust.sh RegCovPheMixedAncestry_split`. This must be done even with a single trait type due to file naming. (sorry)
    - This outputs two phenotype files `Regenie/CovPheno/phenotypesBT.txt` and `Regenie/CovPheno/phenotypesQT.txt` and makes new snplist files that are mac filtered at 1 for the cases that are included in this step to use for regenie step 1.
    - TO CHECK: I have no idea if this step runs if your sample list only has binary or quantitiate traits.
    Now we get into the meat of running regenie
5. Regenie [Step 1](https://rgcgithub.github.io/regenie/overview/#step-1-whole-genome-model): This uses the mac 1 filtered snplists we created to fit a whole genome regression model that captures a good fraction of the phenotype variance attributable to genetic effects.
    - This is done with `./Scripts/collapsing_lclust.sh RegStep1MixedAncestryBT ` or `RegStep1MixedAncestryQT`
    - Sometimes you may run into issues with SNPs being too rare and leading to errors, but these can be removed by creating files called: `Results/Regenie/allQCBAFilter/BT.exclude` and/or `Results/Regenie/allQCBAFilter/QT.exclude` that have the SNP name on each line and including the flag in the regenie call in `Scripts/collapsing_lclust.sh` section: `RegStep1MixedAncestryBT`/`RegStep1MixedAncestryQT` `--exclude $PROJECT/Regenie/allQCBAFilter/BT.exclude \`
6. From here we can switch to dev3-5 for the remaineder of the analysis
    From here we can do actual analyses using regenie.
7. [Step 2 Single variant association testing](https://rgcgithub.github.io/regenie/overview/#step-2-single-variant-association-testing)
    - This is done with `./Scripts/collapsing_lclust.sh RegStep2ExWASBT ` or `RegStep2ExWASQT`
    - Again there is a lot of customization available that can be accessed by changing flags in the `Scripts/collapsing_lclust.sh` file. Please see the [regenie documentation](https://rgcgithub.github.io/regenie/options/#basic-options) for more details.
    - The default setup is an approximate fast firth for binary traits, and LOOCV to deal with case:control imbalance.
    - This uses the minMac setting from the yaml file.
    - The output of regenie files are `*.regenie` and are actually tsv files that can be opened in excel or imported into R with `data.table::fread()`
    - I do not have sophisticated analyses built for the single variant association tests but manhattan plots and qqplots can be made using the qqman library in R, or your library du jour.
    - You can also run saddle point approximation for effect and P-values instead of the fast first using `./Scripts/collapsing_lclust.sh RegStep2ExWASBT_SPA` and `RegStep2ExWASQT_SPA`
8. [Step 2 Gene-based Testing](https://rgcgithub.github.io/regenie/overview/#step-2-gene-based-testing)
    - This is done with `./Scripts/collapsing_lclust.sh RegStep2VCtestsBT` or `RegStep2VCtestsQT`
    - This uses the minMac and aafs setting from the yaml file.
    - Again there is a lot of customization available that can be accessed by changing flags in the `Scripts/collapsing_lclust.sh` file. Please see the [regenie documentation](https://rgcgithub.github.io/regenie/options/#gene-based-testing) for more details.
    - The default setup is an approximate fast firth for binary traits, and LOOCV to deal with case:control imbalance.
    - Default runs the 3 different mask creating schemes (compound het(recessive/comp het model), sum(additive, no max), max(Dominant model))
    - Default runs: `skato-acat,acato-full` [models](https://rgcgithub.github.io/regenie/options/#skatacat-tests) which inlcudes the standard burden `MODEL = ADD-LR` and anumber of ACAT and SKAT tests.
    - This then runs an R script that breaks each master file down by mask (variant type and max minor allele frequence) and creates qq and manhattan plots for each.
    - The default aafs used are: `0.1,0.05,0.01,singleton,all` and the default masks are `LOF, Missense, All, MTR50, MultiPath, and LIMBR50`
    - Of note the output of the `sum` mask type outputs a column called `LOG10P` that has been converted to P-values in the mask subsetted outputs through: `Pval = 10^(-(sum$LOG10P))`
    - - You can also run saddle point approximation for effect and P-values instead of the fast first using `./Scripts/collapsing_lclust.sh RegStep2VCtestsBT_SPA` and `RegStep2VCtestsQT_SPA`
9. Gene-based testing Split.
    - This is done with `./Scripts/collapsing_lclust.sh RegStep2VCtests_QQManhattan`
    - This will split the full regenie files into samller regenie files per phenotype, mask, and QV counting scheme, and will create QQ plots and manhattan plots for each of the same.
    - For SPA results can run the same with: `./Scripts/collapsing_lclust.sh RegStep2VCtests_QQManhattan_SPA`
    
Joshs modified regenie. right now, only works with binary traits. This assumes that you have created a combined ancestry cluster and run coverage on it. Please run the following steps...

clusteringMixedAncestry

clusterCoverageMixedAncestry

RegListVCFMixedAncestry

RegPLINKMergeMixedAncestry_1

RegPLINKMergeMixedAncestry_2

RegCovPheMixedAncestry_1

RegListVCFMixedAncestry_genotype

RegCovPheMixedAncestry_2

RegCovPheMixedAncestry_3

RegCovPheMixedAncestry_split

RegListVCFMixedAncestry_create_AAF

RegStep1MixedAncestryBT

RegStep2VCtestsBT

RegStep2VCtests_QQManhattan

