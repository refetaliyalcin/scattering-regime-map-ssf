# scattering-regime-map-ssf
MATLAB code that draws scattering regime map using static structure factor. Scattering regime map shows the region where independent scattering assumption is valid.

code runs by running main.m

n_pigment = real refractive index of pigment
k_pigment = imaginary refractive index of pigment
n_medium = real refractive index of medium
nang = resolution of scattering angles, should be high for high size parameters, used for numerical convergence, important for calculation of g_dependent.
no_of_f_v = resolution of volume fraction to investigate between 10^-4 and 0.75
no_of_xs = resolution of size parameter to investigate between 10^-2 and 10^3 
