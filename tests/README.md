## Test directory
This is a collection of tests for fluxo 

### Test different builds
For testing different builds (options in cmake), do:
```
python build_tests.py -p 2 -withmpi=1 -buildhdf5=1
```
where all options are explained when typing:
```
python build_tests.py --help
```


### Test the code
When testing the first time, build the meshes:
```
cd meshes
```
see [meshes/README.md](meshes/README.md) for details.

For the freestream test, go to the folder
```
cd freestream
```
see [freestream/README.md](freestream/README.md) for details.
```
cd convergence
```
see [convergence/README.md](convergence/README.md) for details.

each run writes a summary into the test folder.


