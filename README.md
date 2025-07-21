# SWM using the combination of PHS + Poly and the cubed sphere projection  
Author Koji Ogasawara   
ver1. 2025/07/10   

## Index
1. Introduction
2. Required software
3. Contents
4. How to parform Experiments
5. Troubleshooting

# 1. Introduction 
　This repository consists of the spherical advection and shallow water equation models using the combination of the cubed sphere projection method and Radial Basis Function (RBF) generated Finite Difference (RBF-FD) method (Ogasawara and Enomoto submitted). It can perform standard experiments (Williamson et al. 1992) cases 1, 2, and 6. Note that discrepancies in errors between Fortran and Python implementations may arise due to different optimization options.  


# 2. Required software  
### 1 )Python  
    python3
    numpy
    scipy
    matplotlib
    numba
    tqdm

### 2) Fortran  
    PETSc


# 4.How to perform experiments
The parameter settings are described in Ogasawara and Enomoto 2025 submitted.  

### Nodes  
　This repository assumes the use of Maximum Determinant (MD) nodes and icosahedral (icos) grids. If using MD nodes, you need to download the data published on the following website: https://web.maths.unsw.edu.au/~rsw/Sphere/MaxDet/.  The downloaded node data must be rearranged using node_ordering.py, which generates knn_md_N*.txt files. Icosahedral grids, created using spherepts (https://github.com/gradywright/spherepts/), are stored in the nodes directory. The qubarue weights of MD are stored in the downloaded data. For icos grids, the weights were created by Reger and Fornberg (2016) (https://jp.mathworks.com/matlabcentral/fileexchange/51214-spherical_quadrature_rbf-quadrature_nodes) and stored in nodes directory.

### a )Python  
1. Configure params.py.  
2. Execute the model adv.py or swev2.py.

### b )Fortran
1. Configure params.py.
2. Execute save_dmatrix.py and interaction_index.py.
3. Modify md_params.F90 to match the settings in params.py.
4. Compile the code using Make to create an executable file, and then run the experiment.

# 5. Troubleshooting  
1. In the partition_xyz module of swev2.py, save_dmatrix.py, and interaction_index.py, a division-by-zero error may occur.  
For readability, vectors containing node data outside the domain are used in the calculations. NaN (Not a Number) is stored in indices not used in the calculations, so there is no impact on the results.

2. Error of PetscCallA()   
    Please rewrite PetscCallA(petscsubroutine()) as call petscsubroutine()  

# License
  This project is licensed under the MIT License, see the license.txt file for details.

# References  
- K. Ogasawara and T. Enomoto: Application of cubed sphere projection to a meshless shallow water equation model. J. Comput. Phys., submited.  
- Jonah A. Reeger and Bengt Fornberg, 2015: Numerical Quadrature over the Surface of a Sphere. Stud. Appl. Math. 137, 174–188. 
- R.S. Womersley, I.H. Sloan, Interpolation and cubature on the sphere [website], http://web.maths.unsw.edu.au/ rsw/Sphere/ , 2003.   
- Williamson et al. 1992: A Standard Test Set for Numerical Approximations to the Shallow Water Equations in Spherical Geometry,  J. Comput. Phys.,  102, 211-224.  
- Grady. B. Wright, SpherePts [software], https://github.com/gradywright/spherepts/, 2018.
