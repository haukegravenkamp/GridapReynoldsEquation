# GridapReynoldsEquation

**A finite-element solution of the Reynolds equation with a cavitation model.**

This module implements the stabilized finite-element approach published in [1] in the finite-element framework *Gridap* [2] written in Julia. The code complements the paper and can be used to reproduce the results presented there. 
Whenever possible, it uses the same notation and refers to the equation numbers so that the code and paper can be easily followed.

Note: The results published in the paper were not created by this code but by a prototype in Matlab. Hence, **small** differences in the computed errors and residuals can occur when compared with the figures in the paper. 
Also, the shock-capturing approach described in the paper has not yet been included in the Gridap implementation.

I closely followed the Gridap tutorials at [github.com/gridap/Tutorials](https://github.com/gridap/Tutorials) to implement this application. They are great if you want to get started using Gridap.

## content
The module consists of only one file: 
- 'GridapReynoldsEquation.jl' contains the definition of the partial differential equation, domain, finite-element spaces, as well as functions for plotting and performing a convergence test.
  
The other files are related to the specific numerical examples as presented in the paper:

- 'example*.jl' are the input files and main scripts to run in order to reproduce the examples
- 'manufacturedSolutionsMatlab' contains the functions defining body loads. They have been computed using the Matlab program *Manufactured solutions* [3]
- 'results' contains the 'state' files for plotting the results in paraview (optional).

## prerequisites
- Julia
- Paraview (recommended) for plotting the results

## usage
**running the examples:**

Provided you have Julia installed, you should be able to simply download this repository and run the example files.

Details for those new to Julia:
- download or clone this repository
- open a terminal in the main folder GridapReynoldsEquation and start Julia
- start the package manager by typing: ]
- type: instantiate
- wait while Julia installs all missing packages (may take a while)
- run one of the example* files

Running the code for the first time takes significantly longer due to compilation, which is normal behavior in Julia. The second time takes about 20 seconds on my desktop computer to run the complete convergence test for examples 1 or 2.  

**plotting results**
The meshes and results are written in .vtu format in the subfolders of the 'results' folder. They can easily be plotted using Paraview. If you like, you can use the "state" files (.pvsm) provided in the results folder. 
Load them using File->load state, choose 'Search files under specified directory', and make sure to select the directory containing the corresponding .vtu files.

## references  
[1] Gravenkamp, H., Pfeil, S., Codina, R., Stabilized finite elements for the solution of the Reynolds equation considering cavitation, Computer Methods in Applied Mechanics and Engineering (2023)

[2] Gridap, https://gridap.github.io/Gridap.jl/stable/

[3] Gravenkamp, H., Manufactured solutions (2023). doi:10.5281/zenodo.8315789.
20
