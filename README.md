# Umbrella sampling with LAMMPS 
An attempt to use umbrella sampling with LAMMPS.
I used the [LAMMPS LiveCoMS 2025 tutorial](https://github.com/lammpstutorials/lammpstutorials-article) and the [reference paper](doi.org/10.33011/livecoms.6.1.3037).

The workflow also various dependencies. These can be installed using the top-level `environment.yml` file using `micromamba`, `conda` etc.
In addition, WHAM needs to be installed as well; see details below.

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

## Handling dependencies with `micromamba`

The workflow is written in [`Snakemake`](https://snakemake.readthedocs.io/en/stable/) and has certain dependencies. These can be installed in an environment in the following manner: 

```bash
micromamba create -f environment.yml # the first time 
micromamba activate fittingenv # every time you want to activate the environment
```

The environment should then have all the required dependencies, except the WHAM code, which must be installed as described above. 

## Snakemake Workflow 

In a departure from the original tutorial, a separate simulation is run for each umbrella window. In order to accomplish this, for each window, an equilibrated configuration is created first, where the initial position of the particle being pulled is set to the umbrella center. 

Inside the `umbrella_sampling_workflow`, the tree looks like so: 

umbrella_sampling_workflow
├── example_lj
│   ├── config
│        ├── general_config.yml
│        ├── test_samples.csv
│   ├── generate_samples.py
│   ├── results
│   └── runme.sh
└── workflow
    ├── scripts
    └── snakefile

Parameters for the workflow can be changed within the config file (called `general_config.yml`). These are user-defined parameters for the workflow, defining how long to run each production run, umbrella spring force constant value, etc. 

Each window, and all per-window parameters are defined in `test_samples.csv`. This can be generated using `generate_samples.py` and/or modified as well for more flexibility. 

In order to perform a dry run of the workflow,

```bash
cd umbrella_sampling_workflow/example_lj
snakemake --cores 1 --config mpi_threads=1 --snakefile ../workflow/snakefile --directory . -n
```

If the DAG looks OK, you can run the 

