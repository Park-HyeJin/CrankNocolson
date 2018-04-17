By changing the function "Dd_fun", various density-dependent diffusion function can be examined.
Current code implements g(f, rho) = (1-rho)*rho*f*(1-f), bacterial diffusion in the main text. 
If you want to use other diffusion functional form, just change the density dependent part in the function.

To use the program, rd, one needs two argumetns: multiplication factor r and ratio of diffusion coefficient D. 
If there is a configuration at time t in data folder, we can use this configuration usign third argument.
All configuration every 100 time are saved in data folder at a gvien r and D until t=10000 (controled by variable "tmax" in the rd.cpp).
Note that you needs to create "data" folder in advance to save the data. 

Note that be careful to choose your dt. 
In the code we choose dt=0.005.
However, dt has to be larger than (dx)^2/(sigma+1)/4 according to Von Neumann stability analysis, where the factor 4 comes from 2*dimension.
Hence, when you use large sigma value, you need more smaller dt.
