%% 
clc
clear
% close all
%% Define given variables
c_0 = 344; %m/s
rho_0 = 1.2; %kg/m^3

R = 100; %m
f = [3 100 3000 10000]; %Hz

a = 0.1; %m
S = 4*a^2; %m^2

v0 = 1;%m/s
v = [v0, 2*v0, v0, -v0, -2*v0, -v0]; %Velocities %m/s
r_s = [-5*a, -3*a, -a, a, 3*a, 5*a];
q = zeros(size(v)); %volume flow holder

omega = 2* pi*f;
k = omega./c_0;
%% Mesh
x_mesh = -R:R;
y_mesh = 0:R;
[X,Y] = meshgrid(x_mesh,y_mesh);
p_mesh = zeros(length(omega),length(v),length(x_mesh),length(y_mesh));
%% Geometry
phi = 0:0.001:pi;
% x = R*cos(phi);
% y = R*sin(phi);
q = v*S;

for j=1:length(omega)
for i=1:length(v)  
for x = 1:length(x_mesh)
for y = 1:length(y_mesh)
%p(j,i,:) = Piston_Pressure(rho_0,x,y,r_s(i),omega(j),v(i),S);
p_mesh(j,i,x,y) = Piston_Pressure(rho_0,x_mesh(x),y_mesh(y),r_s(i),omega(j),q(i),k(j));
end
end
end
end

%% Plotting setup
phi = -pi/2:0.001:pi/2;
R =100;
yLabel = 'Y Direction m';
xLabel = 'X Direction m';
%%
N_cases = 6;
p_cases = zeros(N_cases,length(omega),length(x_mesh),length(y_mesh));
figName = strings(length(omega),1);
figs = gobjects(length(omega),1);
for i = 1:N_cases
p_cases(i,:,:,:) = sum(p_mesh(:,1:i,:,:),2);
Lp_cases = squeeze(20*log10(abs(p_cases)));
end
for i = 1:length(f)
    figName(i) = strcat('f=', num2str(f(i)),' Hz');
    figs(i) = figure(Name=figName(i),Position =  [100, 0, 1080, 780]);
    for j = 1:length(v)
        subplot(3,2,j)
        contourf(X,Y,squeeze(Lp_cases(j,i,:,:))')
        hold on
        plot(R*sin(phi),R*cos(phi),'r',LineWidth=1.5)
%plot(phi/pi,squeeze(Lp_cases(j,i,:)));
xlabel(xLabel)
ylabel(yLabel)
cb = colorbar;
ylabel(cb,'L_p dB','Rotation',90);
    end
%Save Plots
saveFolder = fullfile(pwd,'\Plots\');

    fileName = strcat('C3_Mesh,',figName(i),'.png');
    filePath = fullfile(saveFolder, fileName);
    exportgraphics(figs(i),filePath,"ContentType","image",'Resolution',600);
end




