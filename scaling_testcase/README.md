## Build the Meshes

A series of meshes are already part of this folder.

If you want to generate others, you should have the most current 
hopr (https://github.com/fhindenlang/hopr.git) installed to build theses meshes!
**To build other meshes**, use the buildmesh python scripts as: 
```
python build_conform_meshes.py <pathTo-HOPR/build/bin/hopr> preproc_hopr.ini
```

## Testcase

Fluxo should be compiled with ```linearscalaradvection``` equation (during cmake), and PARABOLIC=OFF
as the parameterfile is valid only for this equation.

Execute fluxo with the parameterfile as
```
mpirun -np 2 <pathTo-fluxo/build/bin/fluxo> parameter_linadv.ini
```
Parameters than can be changed for the scaling tests are marked in the parameterfile with 
```
! [CAN BE CHANGED]
```
Several meshes are already available in this folder, so to change the mesh, simply change the meshfile parameter:
```
meshfile = CARTBOX_PERIODIC_008_008_008_mesh.h5 
```
The number of elements would be then ```8x8x8=512```.

The maximum number of timesteps can be easily changed with the parameter
```
maxIter = 1000
```
but only if the end time is not reached before, so set tend to a high value, to control the number of timesteps with maxIter! 
```
tend    = 100000.
```
For performance measurements, its recommendable to also switch off the 
analyze time levels (set ```=tend```), which also switches off any IO
```
Analyze_dt  = 100000.
```

Depending on the chosen polynomial degree ```N```, the number of elements in the mesh, and the number of ranks,
the  the wall-clock time changes.

Therefore, additionally to maxIter, we can also set the maximum wall-clock time, in seconds.
```
maxWCT = 60  !maximum wall-clock time in seconds (default=-1, off)
```
This will only be checked after maxIter is reached. If then, the current WCT < maxWCT, maxIter will be changed to reach maxWCT.

**This way, the code runs at least maxIter, and then checks if it runs until maxWCT.**

### Performance index 

The performance index (PID) is always computed by fluxo.

The PID is given at the each analyze level/final time (called "CALCULATION TIME PER TSTEP/DOF ") and is computed as

```
PID= (wall clock time)*nRanks/((N+1)**3*nElems*maxIter)
```

From experience, the PID depends mainly on the architechture, 
and does not vary much between different polynomial degrees.
For theoretically ideal strong scaling,  the PID should be a constant!
From experience, ideal strong scaling is found as long as ```>2000 DOF``` per rank are used.

This translates to a a limiting number of elements per rank, depending on the polynomial degree:

| ```N```    | ```(N+1)^3```        |```>=2000 DOF/rank``` |
|:----------:|---------------------:|---------------------:|
| ```2```    | ```  27 DOF/elem```  | ``` 74 elems/rank``` |
| ```3```    | ```  64 DOF/elem```  | ``` 32 elems/rank``` |
| ```4```    | ``` 125 DOF/elem```  | ``` 16 elems/rank``` |
| ```5```    | ``` 216 DOF/elem```  | ``` 10 elems/rank``` |
| ```6```    | ``` 343 DOF/elem```  | ```  6 elems/rank``` |
| ```7```    | ``` 512 DOF/elem```  | ```  4 elems/rank``` |
| ```8```    | ``` 729 DOF/elem```  | ```  3 elems/rank``` |
| ```9```    | ```1000 DOF/elem```  | ```  2 elems/rank``` |

so that for example for ```N=4``` and a mesh with ```8x8x8``` elements, 
the number of ranks where the strong scaling leaves the ideal behavior should be found for nRanks>32.

If the number of ```elems/rank``` is lower than the limit, 
it seems that communication latency starts reducing the efficiency in the strong scaling.