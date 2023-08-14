## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Coverage of Case Control Cohort
* Run coverage for entire cohort
* This is done so that when we next remove related individuals, we remove poorly covered over well covered

**Must be run on qs1**

Command: `./Scripts/collapsing_lclust.sh initialCoverage`


### Input:
* Automatically detects the sample file created in the [SampleFileCheck](https://github.com/igm-team/ClusteredCollapsing/wiki/SampleFileCheck) step found in /Data/*.ped.txt

### Output:
* /Results/Coverage
  * txt files of ATAV logs
  * csv files of coverage summaries
  * txt files listing existing samples and exons