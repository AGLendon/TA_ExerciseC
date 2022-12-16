%% 
clc
clear
close all
%% Define given variables
c_0 = 344; %m/s
rho_0 = 1.2; %kg/m^3

R = 100; %m
f = [3 100 3000]; %Hz

a = 0.1; %m
S = 4*a^2; %m^2

v0 = 1;%m/s
v = [v0, 2*v0, v0, -v0, -2*v0, -v0]; %Velocities %m/s
r_s = [-5*a, -3*a, -a, a, 3*a, 5*a];
q = zeros(size(v)); %volume flow holder

omega = 2* pi*f;
k = omega./c_0;
%% Geometry
phi = 0:0.001:pi;
x = R*cos(phi);
y = R*sin(phi);
q = v*S;

p = zeros(length(omega),length(v),length(phi));
sum_p = zeros(length(omega),1,length(phi));
p_tot = zeros(length(omega),length(phi));
for j=1:length(omega)
for i=1:length(v)  

p(j,i,:) = Piston_Pressure(rho_0,x,y,r_s(i),omega(j),q(i),k(j));

% sum_p(j,1,:) = sum_p(j,:,:)+p(j,i,:);
% Lp_sum = 20*log10(abs(sum_p));
% figure;
% plot(phi,squeeze(Lp_sum(j,1,:)));

end
end

%% Plotting setup
yLabel = 'L_p dB re 1 Pa';
xLabel = '\phi \pi radians';
%%
N_cases = 6;
p_cases = zeros(N_cases,length(omega),length(phi));
for i = 1:N_cases
p_cases(i,:,:) = sum(p(:,1:i,:),2);
Lp_cases = 20*log10(abs(p_cases));
end
for i = 1:3
    figure;
    for j = 1:6
        subplot(3,2,j)
plot(phi/pi,squeeze(Lp_cases(j,i,:)));
xlabel(xLabel)
ylabel(yLabel)
    end
end
%%


% 
% Lp = 20*log10(p);
% Lp_sum = 20*log10(sum_p);
% figure;
% plot(phi,squeeze(Lp_sum(1,1,:)));

% %% Cases
% N_cases = 6;
% p_case = zeros(length(omega),N_cases,length(phi));
% for i = 1:N_cases
%     p_case(:,i,:) = sum(p(:,1:i,:),2);
%     Lp_case = 20*log10(p_case);
%     figure;
%     plot(phi,squeeze(Lp_case(1,i,:)))
% end