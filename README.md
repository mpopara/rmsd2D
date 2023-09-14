# rmsd2D
Python script for pair-wise RMSD comparison between two ensembles of protein structures


## General description

This script quantifies similarity between two ensembles of conformational models. In the example provided below, compared are two protein ensembles 
obtained from all-atom MD simulations using two distinct force-fields. Comparison between the ensembles is done in terms of RMSD values using the coordinates 
of the backbone C<sub>&alpha;</sub> atoms. Analysis assumes that the two ensembles are clustered, and takes as an input cluster representatives for each of the 
ensembles and their corresponding population fractions (weights). Only the relevant clusters are considered in the analysis by comparing the most populated clusters
in the ensemble, whose number is user-defined through the _cutoff_ parameter (here _cutoff_ = 50). Two-dimensional C<sub>&alpha;</sub>-RMSD matrix is visualized as an image
with side panels which contain the weights of cluster representatives in the descending order. 
In this particular example, analysis revealed that two ensembles are considerably dissimilar, where on average the smallest RMSD across all pairs of structures is ~14 Angström.
Magnitude of this deviation becomes more clear upon mentioning that RMSD threshold for grouping the structures in one cluster is typically up to 4 Angströms.


![RMSD_FF14SB_vs_FF99SB-disp_dpi300](https://github.com/mpopara/rmsd2D/assets/40856779/b6ab868f-285f-4384-9cd7-c34dba3f3707)

Be aware that such comparison of the ensembles may be inapropriate in some scenarios. It has been shown that for ensembles which are optimized against averaged experimental data, 
comparison of individual structures at an atomistic level is an ill-defined problem, since many different solutions can fit the same averaged data.<sup>1</sup>
In such case, it is advised to use more robust model representations for the comparison of the ensembles, such as inter-residue distograms or 3D density maps.<sup>1</sup>

## Input file requirements

As input files required are:

* cluster representatives for the two ensembles, each provided as a separate .dcd trajectory, although any other trajectory format supported by mdtraj will work as well (.xtc, .nc...).
* topology file provided as .pdb file
* .dat file containing weights (population fractions) of ensemble members for each of the ensembles. This space-delimited file is of a size N<sub>conformers</sub> x 2, where the first column contains indices of the ensemble members,
 and the second column contains their corresponding weights. This script assumes that the order of cluster representatives in the trajectory file follows the same order as in the weights file.
It is not required that the cluster representatives and their corresponding weights are sorted in the descencing order according to their weight. Sorting will be done as part of the analysis. 

## Dependencies

RMSD_matrix.py is a python script built on Python 3.8.8. Script was tested under the following configuration:

* Windows 10
* Python 3.8.8
* mdtraj 1.23.0
* numpy 1.23.0
* matplotlib 3.7.1

## References
1. Dittrich, J.; Popara, M.; Kubiak, J.; Dimura, M.; Schepers, B.; Verma, N.; Schmitz,
B.; Dollinger, P.; Kovacic, F.; Jaeger, K. E.; Seidel, C. A. M.; Peulen, T. O.; Gohlke, H.,
Resolution of Maximum Entropy Method-Derived Posterior Conformational Ensembles of a
Flexible System Probed by FRET and Molecular Dynamics Simulations. J Chem Theory
Comput 2023, 19 (8), 2389-2409.


## Authors

* Milana Popara
