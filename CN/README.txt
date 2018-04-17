Here, we implement the two-dimensional Crank-Nicolson algorithm with secon-order Neuman boundary condition (no-external flux).
To employ diffusion and reaction, we use D'yaknov method.
We adopt the theroy for numerical calculations in http://www.cems.uvm.edu/~tlakoba/math337/index.html.


To use the program, rd, one needs two argumetns: multiplication factor r and ratio of diffusion coefficient D. 
If there is a configuration at time t in data folder, we can use this configuration usign third argument.
All configuration every 100 time are saved in data folder at a gvien r and D until t=10000 (controled by variable "tmax" in the rd.cpp).
Note that you needs to create "data" folder in advance to save the data. 


The gnuplot script for plotting time evolution of the configuration is available (the file gnu.TimeEvolution.gpi).
