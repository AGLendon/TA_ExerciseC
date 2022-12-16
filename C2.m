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

omega_c = sqrt(c_0^4*m/B);  %rad/s   critical frequency

f_c=omega_c/(2*pi);
%% b
f = 500;        %Hz         -Driving frequency
omega = 2*pi*f; %rad/s      -Driving frequency

k_0 = omega/c_0;
k_B = (m/B*omega^2)^0.25;

%K_0 > k_b eq 10.42
theta = asind(k_B/k_0);
figure;
rangek0= 0:0.001:k_0;
polarplot(theta+zeros(size(rangek0)),rangek0)
%%

