## running tests
You need to have the current hopr (https://fhinde@bitbucket.org/fhinde/hopr.git) installed to build theses meshes!

Add the hopr executable as a link into the meshes directory
```
cd meshes
ln -s <pathTo-HOPR/build/bin/hopr> hopr
```
Use the buildmesh python script as: 
```
python build_conform_meshes.py hopr parameter_hopr.ini
python build_mortar_meshes.py  hopr parameter_mortar_hopr.ini
python build_mortar_meshes.py  hopr parameter_mortar_periodic_hopr.ini
```

