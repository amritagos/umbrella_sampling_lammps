#!/bin/bash

micromamba activate fittingenv
snakemake --cores 1 --printshellcmds --keep-incomplete --keep-going --config mpi_threads=1 --snakefile ../workflow/snakefile --directory .