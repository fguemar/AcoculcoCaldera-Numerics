Finite Volume implementation to solve the heat transfer equation in the Acoculco Caldera, Puebla, Mexico. The numerical scheme is implemented in FORTRAN 90. The linear problem is solved with a Gauss-Seidel iteration for which the residual is computed in parallel. This code can be compiled and executed in any operative system with a FORTRAN compiler (OpenMP libraries are required for compilation). A multi-core processor and 8 GB Ram memory are recommended. The program requires three input files given in text format: 

1. "in_CurieDepth.txt" 
2. "in_elevAcoculco.txt"
3. "in_tLogs.txt" 

Files 1 and 2 correspond to the estimated Curie isotherm and the Digital Terrain Elevation Model (DTEM) of the study area. These files are read by the program in order to define constant temperature boundary conditions at the top and bottom boundaries of the three-dimensional domain.

Input file 3 contains the temperature logs of the existing exploratory wells in the area, EAC-1 and EAC-2; first column is temperature (Celsius) and second column is depth below surface. Temperature log of well EAC-1 is given form row 1 to 18, the remaining data correspond to well EAC-2.  

The program generates five output files (where '#' represent a depth (of the local heat sources) in meters below sea level):

1. "temp3D_#.txt" 
2. "temp2D_#.txt"
3. "eac1_#.txt"
4. "eac2_#.txt"
5. "deviation_eac1_#.txt"


Files 1 and 2 are the resulting temperature fields of the simulation. The first file is the temperature field for the upper 6000 m of the 3D domain. File two contains a 2D section of the temperature field at x=589937 m. These files can be post-processed for visualization (e.g. the open source programs ParaView and Gnuplot are options for 3D and 2D visualizations respectively). Files 3 and 4 contain the simulated temperature gradients at the location of wells EAC-1 and EAC-2, the gradients are presented with depths in meters below surface. Finally, file 5 contains the overall temperature difference between measured and simulated temperatures for well EAC-1 (this error is not computed for well EAC-2 since the temperature log for this well is a transient temperature profile).


