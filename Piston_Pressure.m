function [p] = Piston_Pressure(rho_0, x, y, r_s, omega, q,k)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    x_r = x-r_s;
    r = sqrt(x_r.^2+y.^2);
    p = zeros(size(r));

    p = 1j*omega*rho_0.*(q./(2*pi.*r)).*exp(-1j*k*r);
end