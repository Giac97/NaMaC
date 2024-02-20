# NaMaC

NaMaC, Nanoporous Material Characterizer is a compact code devised to perform characterization of nanoporous materials, its main functionalities comprise, so far, a code that computes the coordination number and generalized coordination number from an xyz file and from that it identifies the atoms that lie on the surface, then it computes the porosity of the system by using a Monte Carlo integration method. 

## Surface Identification

The surface identification is done by computing first the coordination number for each atom in the system, beside the coordination number a list of the indeces of the atom's first 12 neighbours within a give cutoff is generated. This list is then used to compute the generalized coordination number as:
$$
GCN_i = \sum_{j \in neigh} \frac{CN_j}{CN_{max}},
$$
where $CN_{max}$ is the maximum value of the coordination number, 12 for fcc structures. The surface atoms are then identified as those having a GCN below a certain threshold value.


