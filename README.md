# DesyTauAnalysesRun2_25ns

Instructions and tips for newcomers

The master branch is used for most analyses while other branches are used for NTuple production or other operations which are heavily dependent on the CMSSW version.

To work with this repository please install the appropriate CMSSW version and follow the instructions in our working twiki:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/DesyTauAnalysesRun2

To upload your code:
- please fork the master branch (or the appropriate branch if needed) of this repository, by clicking the fork button in the top right corner (if working via browser)
- from terminal go in the directory `${CMSSW_BASE}/src/DesyTauAnalyses` and link your personal repository by running `git remote add upstream https://github.com/GITHUB_USERNAME/DesyTauAnalysesRun2_25ns.git` (upstream will be the name of the remote, any other name is valid)
- checkout your working branch with `git checkout from-CMSSW_VERSION` (or check the name of the current working branch with `git branch` in case you do not want to make changes to the master branch)
- stage the changes you want to upload with `git add list-of-edited-files`
- commit your changes adding a small comment with `git commit -m "comment" list-of-edited-files`
- push your changes with `git push upstream from-CMSSW_VERSION`
- from browser make a pull request (PR)
