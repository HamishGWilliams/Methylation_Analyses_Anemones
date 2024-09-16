# Methylation_Analyses_Anemones

## Background

Created a new GitHub repository to help consolidate the methylation analyses scripts and some of the results to make it easier to track between versions and to be able to cloud store the scripts and results together into a central place which can then have access to both the HPC cluster and R, making it easier to transfer and store data and information needed which is shared between these different softwares

I have also downloaded the community version of PyCharm to act as a code manager. I am currently learning how to use python, which uses PyCharm as a code manager in the course, therefore it seems appropriate to then include this particular software as a code manager, which means that I can easily look at and edit my code without needing to have direct access to the HPC cluster on my personal laptop to do so!

## Important Note

There will be several folders in this repository. These will be for the 'script' and 'results'. I will not be storing the data files on this repository as they are too large for either the main GitHub (100MB max) or standard Git Large File Storage (Git LFS - max of 1GB for free accounts). Therefore all data files and non-final results files will be stored on the HPC (maxwell). Figures and other summary statistic data files will be moved manually to the github repository, which can then be accessed freely. 

## Commands to remember
### HPC Cluster

**REMEMBER TO BE IN THE REPOSITORY TO MAKE CHANGES TO IT! ALWAYS PULL BEFORE PUSHING!!**

- git pull/push main origin (pulls/pushes from/to main branch)
- git add . (adds all files to the commit)
- git commit -m "message" (adds a message to the commit)
