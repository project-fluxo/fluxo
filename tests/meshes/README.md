## Build the Meshes
You need to have the most current hopr (https://github.com/fhindenlang/hopr.git) installed to build theses meshes!

| Parameterfile                             | Element with mapping      | Boundary Conditions | Description                             :|
|:-----------------------------------------:|:-------------------------:|:-------------------:|:----------------------------------------:|
|`parameter_CartBoxPeriodic.ini`            | Cartesian elements        | fully periodic      | domain is a box [-1,1]^3     |
|`parameter_ConformBoxTrilinear.ini`        | trilinear element mapping | fully periodic      | domain is a deformed box [-1,1]^3 |
|`parameter_ConformBoxCurved.ini`           | curved mapping (degree 3) | fully periodic      | domain is a deformed box [-1,1]^3 |
|`parameter_DeformedBoxMortar.ini`          | curved mapping (degree from script) | Dirichlet | domain is a deformed box [-1,1]^3, with mortar interfaces |
|`parameter_DeformedBoxMortarPeriodic.ini`  | curved mapping (degree from script) | fully periodic | domain is a deformed box [-1,1]^3, with mortar interfaces |


To build the meshes, add the hopr executable as a link into the meshes directory
```
cd meshes
ln -s <pathTo-HOPR/build/bin/hopr> hopr
```
Use the buildmesh python scripts as: 
```
python build_conform_meshes.py hopr parameter_CartBoxPeriodic.ini
python build_conform_meshes.py hopr parameter_ConformBoxTrilinear.ini
python build_conform_meshes.py hopr parameter_ConformBoxCurved.ini
python build_mortar_meshes.py  hopr parameter_DeformedBoxMortar.ini
python build_mortar_meshes.py  hopr parameter_DeformedBoxMortarPeriodic.ini
```
