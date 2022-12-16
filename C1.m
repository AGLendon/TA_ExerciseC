h = 0.01;   %m
b = 0.02;   %m
S = h*b;    %m^2
rho = 7850; %kg/m^3
l1 = 6.3;   %m
l2 = 5;     %m
E = 2e-11;
B = E*I;

%Sub system n =1:4


% p_d = omega*eta_damp*M(n)*v(n)
c_L = sqrt(E/rho);
c_B = (B/(rho*S)*omega^2);
v = [1;2;3;4] %unknowns

c = [12 13 14; 21 23 24; 31 32 34; 41 42 43]; % rows 1-bending ; 2-long; 3-bendin ; 4-long




P = c.*v.^2