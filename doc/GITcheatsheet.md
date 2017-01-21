#GIT cheatsheet
to clone the code, just type
```
 git clone https://github.com/project-fluxo/fluxo.git fluxo
```
this will copy the repository to the directory fluxo. Then doing a 
```
 git status
```
should tell you something like
On branch master
Your branch is up-to-date with 'origin/master'.

##How to commit local changes:
First check if you are up-to-date with the repository
```
 git fetch 
 git status 
```
If status gives you
```
Your branch is up-to-date with 'origin ...
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
You don't need to push the commit right away to the repository, you can accumulate some commits and then push them all with
```
 git push origin branchname
```
where the branchname should the one you are currently working on, use `git status` to see its name.





##Situations:
 
###A) I want to do a git pull but I already modified locally without adding or committing anything...

You can do a 
```
 git stash
```
this will save the local modifications, then pull and then add the modifications with
```
 git stash apply
```

###B) I want to do a git pull but I already added files locally...
You should commit your changes and then continue with **C)** 

###C) I want to do a git pull but I already committed files locally...
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
 git pull --rebase
```


##GIT supporting plugins for the command line
One addition to your bashrc is the command completion (if needed), found for example [here:git-completion.bash](https://github.com/git/git/blob/master/contrib/completion/git-completion.bash)
and you copy this to your home folder and add your `.bashrc` the line
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





