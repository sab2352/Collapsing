# Clustered Collapsing

## Collapsing Analysis Workflow

This repo is meant to be cloned for each new collapsing project. The workflow is based on previously published IGM collapsing analyses performed for heart failure (Poyvsil et al. 2021) and common epilepsies (Epi25 Consortium 2021). This was forked from Gundula's original repository and will contain modifications by the current crew of IGM developers.

Developers should use the "dev" branch while non-developers should use the "main" branch. "main" reflects the most recent release although this repo is under active development.

The wiki tab has been removed from GitHub because of IGM-Team permissions. All wiki files can be found within the <a href = "/Tutorial"> Tutorial </a> folder or linked below to each step of collapsing analysis.


Prerequisites

- https://github.com/0k/shyaml (see https://redmine.igm.cumc.columbia.edu/issues/7936 for installation instructions)

## Tutorials and Documentation

## Begin Here - Setup
<ol>
<li> <a href = "/Tutorial/Repository-Setup.md" target = "_blank">Setting up Repository:</a> Let's begin by setting up our working environment for collapsing! This links to instructions on cloning and downloading the GitHub repo. Best located on /projects when possible. Each cohort will require its own copy of the repo.</li>

<li> <a href = "/Tutorial/Yaml-Setup.md" target = "_blank">YAML file:</a> The code uses a configuration file (.yaml). This is the only part of the code that is edited by the user and is how we set up our variables for each step of collapsing analysis. Once done, you're ready to begin!</li>
</ol>

## Analysis 
<ol>
<li> <a href = "/Tutorial/Creating-Cohort.md" target = "_blank">cohortSelection:</a> The first step of collapsing analysis is when we select our samples for our cases and controls. In this step, you create your ped file. This is a file of cases and controls.</li>

<li> <a href = "/Tutorial/SampleFileCheck.md" target = "_blank">SampleFileCheck:</a> Once the cohort is created, this quick step is to check that the samples in our cases and controls are consistent with the samples in the IGM databases and that each row is formatted correctly for ATAV.</li>

<li> <a href = "/Tutorial/Coverage-of-CaseControl-Cohort.md" target = "_blank">initialCoverage:</a> Next, coverage is performed on the case/control cohort. This step of collapsing removes samples that are not adequately covered between cases and controls. The output for this step is used in the Kinship step when removing related individuals (given a choice, the program keep the better covered individual).</li>


<li> <a href = "/Tutorial/Creation-of-Principle-Components-and-Removal-of-Related-Samples.md" target = "_blank">CohortPCA:</a> This step is the Kinship and FlashPCA step. In this step, we create ancestry-based principle components and remove related individuals. This step uses the coverage data from initialCoverage to keep the best covered samples.</li>

<li>Clustering Step</li>
	
<ul>
<li> <a href = "/Tutorial/Create-Clusters.md" target = "_blank">clustering:</a> This step uses the PCAs generated in the prior step and Louvain clustering to create geographic ancestry-based case/control clusters.</li>

<li> <a href = "/Tutorial/Create-Mixed-Ancestry-Cluster.md" target = "_blank">clusteringMixedAncestry:</a> For those interested in analyzing a single cluster with all ancestries, this step facilitates that. It can be used in Regenie analyses.</li>


</ul>

<br>
<h4> After the clustering step, the following are run within the clusters: </h4>
<br>


<li> <a href = "/Tutorial/pcaCluster.md" target = "_blank">pcaCluster:</a> Re-runs ancestry PCA on each cluster. In doing so, excludes ancestry outliers. This step can be skipped if desired.</li>

<li> <a href = "/Tutorial/Coverage-by-Cluster.md" target = "_blank">clusterCoverage:</a> Coverage harmonization is performed within each cluster. This step removes variants that are not adequately covered between cases and controls within each cluster. This step can be skipped if desired.</li>

<li> <a href = "/Tutorial/Master-Variant-File.md" target = "_blank">createMasterGeno:</a> This step creates large variant files (master genotypes files) for each cluster. Subsequent models will be drawn from these files so only one trip to the database is necessary.</li>

<li> <a href = "/Tutorial/Create-Models.md" target = "_blank">models:</a> It's time to select the models that you want to run for your analysis! This step filters the variants in the master genotype files created in the preceding step in order to create specific models. You can read more about the options for each models (or even create your own!) here: <a href = "/Tutorial/selecting_models_options.md" target = "_blank">ATAV options for models</a>. </li>

<br>
<h4> The following steps aggregate the analysis across the included clusters: </h4>
<br>

<li> <a href = "/Tutorial/Gene-Based-CMH.md" target = "_blank">CMH:</a> The Cochran-Mantel-Haenszel (CMH) test allows us to make sense of the significance values of each ancestry cluster together.</li>

<li> <a href = "/Tutorial/CMH-Permutations.md" target = "_blank">CMH Permutations:</a> This script creates the permutations needed for CMH test. The number of permutations can be adjusted in the YAML file.</li>

<li> <a href = "/Tutorial/Create-QQ-Plots.md" target = "_blank">QQ Plot:</a> A QQ plot allows us to check for inflation/deflation and find top hits!</li>

</ol>


## Special Functionalities 
<ol>
<li> Mixed Ancestry Cluster 
The following steps can be performed for the mixed ancestry cluster:

<ul>
<li>pcaMixedAncestry</li>
<li>clusterCoverageMixedAncestry</li>
<li>createMasterGenoMixedAncestry</li>
<li>modelsMixedAncestry</li>
</ul>

</li>

<li> Running Coverage by Sex
 
<ul>
<li>clusterCoverageMale</li>
<li>clusterCoverageFemale</li>
</ul>

</li>

</ol>

## Supplemental Features 
The following are optional steps usually necessary for methods or supplemental data for collapsing papers and presentations. 

<ul>
	<li> <a href = "/Tutorial/Create-Cluster-Map.md" target = "_blank">Create Cluster Map:</a> This allows us to better visualize the clusters that we create.</li>
	<li> <a href = "/Tutorial/Track-Sample-Removal.md" target = "_blank">Track Sample Removal:</a> This allows us to keep track of the number of samples removed in each step.</li>
</ul>

## Other Analyses 
But wait! If you are not interested in traditional collapsing analysis, our repo can also help you run these additional analyses — please take a look! 
<ol>
	<li> <a href = "/Tutorial/Run-Regenie.md" target = "_blank">Regenie:</a> Run Regenie 3.1.3 on your data.</li>
	<li> <a href = "/Tutorial/Run-NHC.md" target = "_blank">NHC</a></li>
	<li> <a href = "/Tutorial/Forest-Plot.md" target = "_blank">Forest Plot:</a> Generate a forest plot if provided a csv with ORs, CIs, P values, etc.</li>
	<li> <a href = "/Tutorial/findThreshold.md" target = "_blank">findThreshold:</a> Detects optimal threshold for gene set analysis. Currently only working for LOEUF.</li>
	<li> <a href = "/Tutorial/create_gene_set_model.md" target = "_blank">Gene Set Model</a></li>
	<li> <a href = "/Tutorial/create_digenic.md" target = "_blank">Create Digenic Model:</a> Create a digenic model out of a single model.</li>
</ol>

## ATAV Wiki
For granular detail regarding ATAV functions, please refer to the <a href = "https://redmine.igm.cumc.columbia.edu/projects/atav/wiki">ATAV wiki</a>.

## FAQ and Troubleshooting 
<ul>
<li> <a href = "/Tutorial/Frequently-Asked-Questions.md" target = "_blank">Frequently Asked Questions</a> </li>
	
<li> <a href = "/Tutorial/Troubleshooting.md" target = "_blank">Troubleshooting</a> </li>

</ul>