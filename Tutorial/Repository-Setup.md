
## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Repository Setup 

1. Create a folder to house your project. Each project is a set of samples from beginning to end. 
2. Save repo in the folder you just created: 
	- Method 1: Easiest way to do this is to go to https://github.com/igm-team/ClusteredCollapsing and click the green "Code" button, download ZIP into your projects folder and unzip the folder.
	
	- Method 2: Another way is to clone the repo. This will allow you to pull updates as the repo gets updated. 

	An example code is below:
	```
	@dev2:/nfs/projects/YOURPROJECT$ git clone https://github.com/igm-team/ClusteredCollapsing
	```

3. Create a .txt file with a list of samples to be included. No header. You will modify the YAML file with a path to this file to create your cohort. 
4. The main bash script (collapsing_lclust.sh) has to be turned into an executable form.
	- You can do this by doing: 
	```
	YOURPROJECT$ cd ./ClusteredCollapsing
	ClusteredCollapsing$ chmod +x ./Scripts/collapsing_lclust.sh 
	```
5. Move collapsing.yaml from the Examples folder to Input folder and then edit the Yaml file with the variables you choose for each step.
	- See [Yaml Setup](/Tutorial/Yaml-Setup.md) for more info

NOTES:
- We recommend projects be based in /nfs/projects to facilitate sharing results and handing off projects from one person to the next.
- A folder on /nfs/projects can be requested on https://redmine.igm.cumc.columbia.edu/projects/storage-requests with the ticket assigned to W. John Burns.
- All commands should be run from within the ClusteredCollapsing folder.
- Rscript commands must be run either on your local machine or on dev1-dev5 while ATAV commands must be run on qs1.
- ATAV output will have a timestamp before the filenames. The collapsing_lclust.sh ignores the timestamp but will get confused if there are more two of the same file but with different timestamps. If you accidentally run a step twice (it happens!), delete the files that you don't need.

