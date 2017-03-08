## running tests
after having built the meshes (in ../meshes) before executing the tests!

If you have compiled with linearscalaradvection equation, use (-p option is number of procs):
```
python runtest.py -p 1 ../../build/bin/fluxo parameter_convergence_linadv.ini

python runtest.py -p 1 ../../build/bin/fluxo parameter_convergence_linadvdiff.ini
```

