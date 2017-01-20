#GIT cheatsheet
to clone the code, just type
```
 git clone https://github.com/project-fluxo/fluxo.git fluxo
```
this will copy the repository to the directory fluxo.

##Situations:
 
**I want to do a git pull but I already modified locally without adding or committing anything...**

You can do a 
```
 git stash
```
this will save the local modifications, then pull and then add the modifications with
```
 git stash apply
```

