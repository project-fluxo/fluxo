# GIT cheatsheet
## GIT clone and push
To copy the code to your local machine, so-called `clone`, just type
```
 git clone https://github.com/project-fluxo/fluxo.git fluxo
```
this will copy the full git repository to the directory fluxo. Change to the fluxo directory.

Doing a 
```
 git status
```
should tell you something like
```
On branch master
Your branch is up-to-date with 'origin/master'.
```

If you want to work on an existing branch instead of the master, just type
```
 git checkout branchname
 git status
```

The list of all local and remote branches is displayed with
```
 git branch -rvva
```

If you want to create a new branch (only locally, must be pushed to add it to the remote)
```
 git checkout -b new-branchname
```
this will copy all files of the branch you were before into the new branch.


## How to pack current code into a tar.gz
This will pack only the files in the repository **excluding** the git repository itself. Execute the command in the main fluxo directory
```
 git archive --format=tar.gz --prefix=fluxo-v1_1/ HEAD > ../fluxo-v1_1.tar.gz  
```
This will put the 'fluxo-v1_1.tar.gz' on the same level as the fluxo directory. 

If you want to pack the code including the git repository, you should exclude 'shared' and 'build' directories...


## How to commit local changes <a name="head_commit"></a>:
First check if you are up-to-date with the repository
```
 git fetch 
 git status 
```
If status gives you
```
Your branch is up-to-date with 'origin'...
```
then you can continue here, else look at the `Situations` below.

If your staging area is empty and you only have local modifications, 
you should choose the files you want to add to the staging area by
```
 git add path-to-file
```
Once you added all files you want to commit (does not have to be all modified files), 
you commit to your local repository with
```
 git commit -m" commit message"
```
Please add a meaning full commit message to your commit.
You do not need to push the commit right away to the repository, you can accumulate some commits and then push them all with
```
 git push origin branchname
```
where the branchname should the one you are currently working on, use `git status` to see its name.


## Situations: "I want to update my local repo with the github remote..."

### A) ... and I have no local changes...
Check if you really have no changes
```
 git status 
```
and update with 
```
 git pull origin currentbranchname
```
 
### B) ... but I already modified locally without adding or committing anything?!?

You can do a 
```
 git stash
```

this will save the local modifications, then pull
```
 git pull origin currentbranchname
```

and then add the modifications with
```
 git stash apply
```

### C) ... but I already added files locally?!?
You should [commit](#head_commit) your changes and then continue with **D)** 

### D) ... but I already committed files locally?!?
Then, first you do a 
```
 git fetch 
```
which basically only updates the information from your repository. doing a 
```
 git status
```
will tell you that the local and the remote have diverged. The common practice is to put the commits of the repository **first**
and add the local changes **on top** by doing a
```
 git pull --rebase
 git status
```
If there is no conflict, you should be ready to push.
If a conflict arises ('both modified'), look into the file and resolve it (marked with '<<<' '>>>' signs). By adding all resolved files, you are ready to finish the process by doing
```
 git rebase --continue
```
and then you should be ready to push.
```
 git push origin branchname
```

## GIT supporting plugins for the command line
One addition to your bashrc is the git command completion (if needed), found for example here [git-completion.bash](https://github.com/git/git/blob/master/contrib/completion/git-completion.bash)
Just copy the file to your home folder and add your `.bashrc` the line
```
source git-completion.bash
```

A more sophisticated tool is the [bash-git-prompt](https://github.com/magicmonty/bash-git-prompt), it helps to keep the overview.
Just copy the directory to your home folder and add to your `.bashrc`
```
source ~/.bash-git-prompt/gitprompt.sh
GIT_PROMPT_ONLY_IN_REPO=1
GIT_PROMPT_THEME=Single_line
```





