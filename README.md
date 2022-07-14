# CuboidSpinfoamScalarField

In this repository, we present the code used to study a scalar field coupled to a restricted 4D spin foam model (a path integral approach to quantum gravity), whose results are contained in the following paper: https://arxiv.org/abs/2206.04076

This algorithm defines a Metropolis algorithm exploring the coupled spin foam - scalar field configuration space. The different ingredients of the code are contained in separate files.

## spin_foam.jl

Contains the definition of semi-classical amplitudes, essentially the weights given to the irregular lattices described by the spin foam model.

## lattice-svec.jl

Definition of the lattice (periodic boundary conditions etc.) and the discrete scalar field action for irregular lattices (derived from disrete exterior calculus).

## proposal-svec.jl

Contains proposal methods for new configurations: on the one hand for lengths of the spin foam, on the other hand scalar field configurations.

## imp_sampling-svec.jl

Contains the actual Metropolis algorithm: thermalization / burn-in and generation of samples. In both cases, we modify a given configuration by randomly deciding whether to modify lengths or scalar field variables. This is accepted or rejected according to how the probability distribution changes compared to the old configuration.

The state at the end of thermalization is stored and can be used as the starting point to begin sampling or do another thermalization run.

## main-svec.jl / main_array-svec.jl

Main file for starting the code. The "array" version works for a range of the parameter alpha (passed as an integer).

In this file, one can change the lattice size (must be greater or equal to 3), thermalization steps, number of samples and steps in between them. Also determines the folder where to store the sampled results for later analysis.

The mass of the scalar field is passed to this function, e.g. in the slurm script.

## plots.jl

Generates the desired plots, from thermalization / sample plots to check whether the system has thermalized, to correlation functions of the scalar fields.