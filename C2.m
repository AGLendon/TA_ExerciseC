%% Q2
clear; clc; close all;
%%
h = 0.05;   %m          -depth of beam
rho = 7850; %kg/m^3     -plate density
E = 2e11;   %N/m^2      -Youngs modulus
mu = 0.3;   %           -Poisson Ratio
c_0 = 344;  %m/s        -speed of sound in air

m = rho*h;
B = (E*h^3)/(12*(1-mu^2));

omega_c = sqrt(c_0^2*m/B);  %rad/s   critical frequency

