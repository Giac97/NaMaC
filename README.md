# NaMaC
<img src="https://github.com/Giac97/NaMaC/assets/126502167/2d04f878-6827-462b-99d9-c61f3d51bffb" alt="namac_logo" width="300" height="200">

NaMaC, Nanoporous Material Characterizer is a compact code devised to perform characterization of nanoporous materials, its main functionalities comprise, so far, a code that computes the coordination number and generalized coordination number from an xyz file and from that it identifies the atoms that lie on the surface, then it computes the porosity of the system by using a Monte Carlo integration method. 

## Surface Identification

The surface identification is done by computing first the coordination number for each atom in the system, beside the coordination number a list of the indeces of the atom's first 12 neighbours within a give cutoff is generated. This list is then used to compute the generalized coordination number as:
$$GCN_i = \sum_{j \in neigh} \frac{CN_j}{CN_{max}},$$
where $CN_{max}$ is the maximum value of the coordination number, 12 for fcc structures. The surface atoms are then identified as those having a GCN below a certain threshold value.

## Porosity Calculation

The code computes the porosity using a simple Monte Carlo integration scheme, the insertion of a probe atom is attempted at a random location within the structure, if the probe is found to overlap with any of the atoms inside the material, the insertion is rejected, else it is accepted. This procedure is carried out for a selected number of samples, the porosity is then found as:
$$\phi = \frac{\text{accepted}}{\text{number of samples}}$$
the probe has a certain radius so that it can stand for an actual atom as those used in experiments such as Helium or Mercury, or have a radius zero to be more of an ideal probe.

## Input File

The input file is a text file containing parameters for the calculation, it is a Fortran namelist with the following structure:
```fortran
&parameters
    fname = "slice.xyz",
    cutoff = 3.5,
    CNmax = 12.0,
    gcn_cutoff = 10.0,
    pbc = 0
    /
&porosity
    n_samples = 100000
    r_probe = 0.0,
    r_atom = 2.5
/
```

The **parameters** inputs are:
- fname: the name of the xyz file containing the structure
- cutoff: the neighbour cutoff distance
- CNmax: The maximum coordination number (12 for fcc)
- gcn_cutoff: The GCN cutoff to separate surface from bulk
- pbc: Whether to use (1) or not (0) periodic boundary codnitions when finding neighbors (to be implemented!)

The **porosity** inputs are:
- n_samples: The numbe rof attempted insertions of probe atoms
- r_probe: The radius of the probe atom
- r_atom: The radius of the atoms composing the system (For now only monometallic systems)

## Compilation

A Makefile is provided with the code, this however requires the PGI Fortran compiler part of the Nvidia HPC toolkit, will be updated to add the possibility to compile for serial CPU machines using gfortran when possible
