% TA assignment 1c 
% question 1
clear;
close all;

E= 2*10^11; % young's modulus of the beam in N/m2
rho= 7850;  % density of the beams in kg/m3
L_h= 6.3;   % length of the horizontal beam in m
L_v= 5;     % length of the vertical beam in m
h= 0.01;    % heigth of the beams in m
b=0.02;     % depth of the beams in m
eta= 0.05;  % damping loss factor
S=b*h;      % cross sectional area of the beam
M_h= rho*S*L_h;  % horizontal beam mass in kg
M_v= rho*S*L_v;  % vertical beam mass  in kg
f= 100:10:3000; % frequency 
omega= 2*pi*f;  % angular freq

% Calculate the speed of longitudinal waves cl= sqrt(E/rho) and group speed cgl=cl;

cl= (E/rho)^(0.5); % speed in long. waves.
cgl= cl;  % group speed 

% v = zeros(1,length(f));
for i = length(f):-1:1
% Calculate the speed of bending waves cb= ((B*omega(i)^2)/(rho*S))^(0.25)
I= 1/12*(b*h^3);
B= E*I;
k_b = ((omega(i).^2*rho*S)/B).^0.25;
cb = omega(i)./k_b;

cgb=2*cb;


%calculate wavelenght for bending waves, lambda_b
lambda_b= cb./f(i);


%calculate tao ( 4 of them are used) (chapter 6)

beta= (2*h)./lambda_b;  %ratio between thicknesses and bending waves' wavelength 
tao_bb= (1+2*beta.^2)./(2+6*beta+9*beta.^2); %bending to bending
tao_bl_lb= (5*beta+8*beta.^2)./(2+6*beta+9*beta.^2); %bending to longititunal or the opposite
rho_bl_lb= beta./(2+6*beta+9*beta.^2); %bending to longitidunal reflected or the opposite
%rho_bb= (1-beta.^2)./(2+6*beta+9*beta.^2); % bending to bending reflected
tao_ll= (beta.^2)./ (2+6*beta+9*beta.^2); % longitidunal to longitidunal
%rho_ll= 2./(2+6*beta+9*beta.^2); %  % longitidunal to longitidunal reflected

%calculate coupling loss factors eta_ij=(tao_ij/c_group_i)/(omega(i)*L_beam*(2-tao_ij))
%h values are identical so beta values is same for the beams.

%for eta_12 use tao_12= rho_bl and cgb; eta_21 use rho_lb and cgl; eta_31 use tao_bb and cgb; eta_41 use tao_lb and cgl
%for eta_13 use tao_13= tao_bb and cgb; eta_23 use tao_lb and cgl; eta_32 use tao_bl and cgb; eta_42 use tao_ll and cgl
%for eta_14 use tao_14= tao_bl and cgb; eta_24 use tao_ll and cgl; eta_34 use rho_bl and cgb; eta_43 use rho_lb and cgl

eta_12= (rho_bl_lb.*cgb)./((omega(i).*L_h).* (2-rho_bl_lb));
eta_13= (tao_bb.*cgb)./((omega(i).*L_h).* (2-tao_bb));
eta_14= (tao_bl_lb.*cgb)./((omega(i).*L_h).* (2-tao_bl_lb));

eta_21= (rho_bl_lb.*cgl)./((omega(i).*L_h).* (2-rho_bl_lb));
eta_23= (tao_bl_lb.*cgl)./((omega(i).*L_h).* (2-tao_bl_lb));
eta_24= (tao_ll.*cgl)./((omega(i).*L_h).* (2-tao_ll));

eta_31= (tao_bb.*cgb)./((omega(i).*L_v).* (2-tao_bb));
eta_32= (tao_bl_lb.*cgb)./((omega(i).*L_v).* (2-tao_bl_lb));
eta_34= (rho_bl_lb.*cgb)./((omega(i).*L_v).* (2-rho_bl_lb));

eta_41= (tao_bl_lb.*cgl)./((omega(i).*L_v).* (2-tao_bl_lb));
eta_42= (tao_ll.*cgl)./((omega(i).*L_v).* (2-tao_ll));
eta_43= (rho_bl_lb.*cgl)./((omega(i).*L_v).* (2-rho_bl_lb));


% Set up the SEA model that includes four subsystem: (M_h=M1=M2; M_v=M3=M4) 

% general formula:
 % P_d1=omega(i)*M1*eta_1*(V_1)^2 
 % P_12= c_12* (V_1)^2 ; P_21= c_21* (V_2)^2 
 % c_12= omega(i)* M1* eta_12

%  subsystem 1: bending waves in horizontal beam B_1
%  P_in=1; P_in + P_21- P_12 + P_31 - P_13 + P_41 -P_14 - P_d1 = 0
% -(omega(i)* M1* (eta+ eta_12 + eta_13 + eta_14)) * (V_1)^2 +  (omega(i)*M2*eta_21*(V_2)^2 + omega(i)* M3*eta_31*(V_3)^2 + omega(i)*M4*eta_41*(V_4)^2 =0

%  subsystem 2: longitudinal waves in horizontal beam B_1
%  P_12 - P_21 + P_32 - P_23 + P_42 - P_24 - P_d2 = 0
%  omega(i)* M1* eta_12* (V_1)^2  -(omega(i)* M2* (eta+ eta_21 + eta_23 + eta_24))*(V_2)^2 + omega(i)* M3*eta_32*(V_3)^2 + omega(i)*M4*eta_42*(V_4)^2=0

%  subsystem 3: bending waves in vertical beam B_2
%  P_13 - P_31 + P_23 - P_32 + P_43 - P_34 - P_d3 = 0
%  omega(i)* M1* eta_13* (V_1)^2 + (omega(i)* M2*eta_23*(V_2)^2 -(omega(i)*M3* (eta+ eta_31 + eta_32 + eta_34))*(V_3)^2 + omega(i)*M4*eta_43*(V_4)^2=0


%  subsystem 4: longitudinal waves in vertical beam B_2
%  P_14 - P_41 + P_24 - P_42 + P_34 - P_43 - P_d4 = 0
%  omega(i)* M1* eta_14* (V_1)^2 + (omega(i)* M2*eta_24*(V_2)^2 + omega(i)*M3* eta_34*(V_3)^2 - (omega(i)* M4* (eta+ eta_41 + eta_42 +eta_43))*(V_4)^2 =0





 A = [ -omega(i)*M_h*(eta+eta_12+eta_13+eta_14),      omega(i)*(M_h*eta_21),                       omega(i)*(M_v*eta_31),                       omega(i)*(M_v*eta_41);
        omega(i)*(M_h*eta_12),                       -omega(i)*M_h*(eta+eta_21+eta_23+eta_24) ,    omega(i)*(M_v*eta_32),                       omega(i)*(M_v*eta_42);  
        omega(i)*(M_h*eta_13),                        omega(i)*(M_h*eta_23),                      -omega(i)*M_v*(eta+eta_31+eta_32+eta_34),     omega(i)*(M_v*eta_43);
        omega(i)*(M_h*eta_14),                        omega(i)*(M_h*eta_24),                       omega(i)*(M_v*eta_34) ,                     -omega(i)*M_v*(eta+eta_41+eta_42+eta_43)];
Y = [-1;0;0;0];
v(:,i) = linsolve(A,Y);
Lv(:,i) = 10*log10(v(:,i)./v(1,i));
end


%%
lambda_b= cb./f;


%calculate tao ( 4 of them are used) (chapter 6)

beta= (2*h)./lambda_b;  %ratio between thicknesses and bending waves' wavelength 
tao_bb= (1+2*beta.^2)./(2+6*beta+9*beta.^2); %bending to bending
tao_bl_lb= (5*beta+8*beta.^2)./(2+6*beta+9*beta.^2); %bending to longititunal
rho_bl_lb= beta./(2+6*beta+9*beta.^2); %bending to longitidunal reflected
rho_bb= (1-beta.^2)./(2+6*beta+9*beta.^2); % bending to bending reflected
tao_ll= (beta.^2)./ (2+6*beta+9*beta.^2); % longitidunal to longitidunal
rho_ll= 2./(2+6*beta+9*beta.^2); %  % longitidunal to longitidunal reflected
%%
f1 = figure(Name= 'Velocity Levels',Position =  [100, 0, 880, 780]);
semilogx(f, Lv); 
legend('L_{v1}','L_{v2}','L_{v3}','L_{v4}')
xlabel('f Hz')
ylabel('Velocity Level dB ref. v_1')
xlim([100 3000])
ylim([-25 1.1])
grid on;
%%

f2 = figure(Name='Transmission & reflection coefficients bend, case:bending',Position =  [100, 0, 880, 780]);
semilogx(beta.^2, tao_bb);
hold on
sum2 = tao_bl_lb+tao_bb;
sum3 = rho_bl_lb+tao_bl_lb+tao_bb;
sum4 = rho_bb+rho_bl_lb+tao_bl_lb+tao_bb;
semilogx(beta.^2, sum2); 
semilogx(beta.^2, sum3); 
semilogx(beta.^2, sum4); 
%legend('tao_{bb}','tao_{bl}','rho_{bl}','rho_{bb}')
xlim([beta(1)^2,beta(end)^2])
ylim([0 1.02])
xlabel('\beta^2')
%ylabel()
grid on;

index = 200;
step = 0.001;
taubbRange = 0:step:tao_bb(index);
vLine = zeros(size(taubbRange));
plot(vLine+beta(index).^2,taubbRange,'k-',HandleVisibility='off')
text(beta(index+5).^2,taubbRange(fix(length(taubbRange)/2)),'\tau_{bb}')

index = 180;
taublRange = tao_bb(index):step:sum2(index);
vLine = zeros(size(taublRange));
plot(vLine+beta(index).^2,taublRange,'k-',HandleVisibility='off')
text(beta(index+5).^2,taublRange(fix(length(taublRange)/2)),'\tau_{bl}')

index = 200;
rhoblRange = sum2(index):step:sum3(index);
vLine = zeros(size(rhoblRange));
plot(vLine+beta(index).^2,rhoblRange,'k-',HandleVisibility='off')
text(beta(index-30).^2,rhoblRange(fix(-10+length(rhoblRange)/2)),'\rho_{bl}')


index = 180;
rhobbRange = sum3(index):step:sum4(index);
vLine = zeros(size(rhobbRange));
plot(vLine+beta(index).^2,rhobbRange,'k-',HandleVisibility='off')
text(beta(index+5).^2,rhobbRange(fix(length(rhobbRange)/2)),'\rho_{bb}')

%%
f3 = figure(Name='Transmission & reflection coefficients Long, case:bending',Position =  [100, 0, 880, 780]);
semilogx(beta.^2, tao_ll);
hold on
sum2 = tao_bl_lb+tao_ll;
sum3 = rho_bl_lb+tao_bl_lb+tao_ll;
sum4 = rho_ll+rho_bl_lb+tao_bl_lb+tao_ll;
semilogx(beta.^2, sum2); 
semilogx(beta.^2, sum3); 
semilogx(beta.^2, sum4); 
%legend('tao_{bb}','tao_{bl}','rho_{bl}','rho_{bb}')
xlim([beta(1)^2,beta(end)^2])
ylim([0 1.02])
xlabel('\beta^2')
%ylabel()
grid on;

index = 200;
step = 0.001;
taullRange = 0:step:tao_ll(index);
vLine = zeros(size(taullRange));
plot(vLine+beta(index).^2,taullRange,'k-',HandleVisibility='off')
text(beta(index+5).^2,0.05+taullRange(fix(length(taullRange)/2)),'\tau_{ll}')

index = 180;
taulbRange = tao_ll(index):step:sum2(index);
vLine = zeros(size(taulbRange));
plot(vLine+beta(index).^2,taulbRange,'k-',HandleVisibility='off')
text(beta(index+5).^2,taulbRange(fix(length(taulbRange)/2)),'\tau_{lb}')

index = 200;
rholbRange = sum2(index):step:sum3(index);
vLine = zeros(size(rholbRange));
plot(vLine+beta(index).^2,rholbRange,'k-',HandleVisibility='off')
text(beta(index-30).^2,rholbRange(fix(-10+length(rholbRange)/2)),'\rho_{lb}')


index = 180;
rhollRange = sum3(index):step:sum4(index);
vLine = zeros(size(rhollRange));
plot(vLine+beta(index).^2,rhollRange,'k-',HandleVisibility='off')
text(beta(index+5).^2,rhollRange(fix(length(rhollRange)/2)),'\rho_{ll}')

%%
f4 = figure (Name='CORRECTED',Position =  [100, 0, 880, 780]);
semilogx(f, v.^0.5); 
xlabel('f Hz')
ylabel('Velcoity ms^{-1}')
legend('v_1','v_2','v_3','v_4')
grid on;

%xlim([beta(1)^2,beta(end)^2])
%ylim([0 1])
thickenall_big

%% Export Figures
saveFolder = fullfile(pwd,'\Plots\');

    fileName = strcat('C3_VelocityLevels','.png');
    filePath = fullfile(saveFolder, fileName);
    exportgraphics(f1,filePath,"ContentType","image",'Resolution',600);

    fileName = strcat('C3_BendingCoefficients','.png');
    filePath = fullfile(saveFolder, fileName);
    exportgraphics(f2,filePath,"ContentType","image",'Resolution',600);

    fileName = strcat('C3_LongitudinalCoefficients','.png');
    filePath = fullfile(saveFolder, fileName);
    exportgraphics(f2,filePath,"ContentType","image",'Resolution',600);

    fileName = strcat('C3_Velocity','.png');
    filePath = fullfile(saveFolder, fileName);
    exportgraphics(f4,filePath,"ContentType","image",'Resolution',600);