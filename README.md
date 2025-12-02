# Umbrella sampling with LAMMPS 
An attempt to use umbrella sampling with LAMMPS.
I used the [LAMMPS LiveCoMS 2025 tutorial](https://github.com/lammpstutorials/lammpstutorials-article) and the [reference paper](doi.org/10.33011/livecoms.6.1.3037).

## Installing WHAM

I used the [WHAM program implemented by Alan Grossman](http://membrane.urmc.rochester.edu/?page_id=126). 

The installation is straightforward, but remember to use $CONDA_PREFIX when installing in the environment. 
```bash
micromamba activate fittingenv
cd wham-directory
mkdir build && cd build
cmake ..
cmake --build .
cmake --install . --prefix $CONDA_PREFIX
```

