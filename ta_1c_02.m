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
M1= rho*S*L_h;
M2= rho*S*L_v;
f= 100:10:3000;
omega= 2*pi*f;

% Calculae the speed of longitudinal waves cl= sqrt(E/rho) and group speed cgl=cl;

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
M = [M1;M2;M1;M2];

%calculate wavelenght for bending waves, lambda_b
lambda_b= cb./f(i);


%calculate tao ( 4 of them are used) (chapter 6)

beta= (2*h)./lambda_b;  %ratio between thicknesses and bending waves' wavelength 
tao_bb= (1+2*beta.^2)./(2+6*beta+9*beta.^2); %bending to bending
tao_bl_lb= (5*beta+8*beta.^2)./(2+6*beta+9*beta.^2); %bending to longititunal
rho_bl_lb= beta./(2+6*beta+9*beta.^2); %bending to longitidunal reflected
%rho_bb= (1-beta.^2)./(2+6*beta+9*beta.^2); % bending to bending reflected
tao_ll= (beta.^2)./ (2+6*beta+9*beta.^2); % longitidunal to longitidunal
%rho_ll= 2./(2+6*beta+9*beta.^2); %  % longitidunal to longitidunal reflected

%calculate coupling loss factors eta_ij=(tao_ij/c_group_i)/(omega(i)*L_beam*(2-tao_ij))
%h values are identical so beta values is same for the beams.

%for eta_12 use tao_12= tao_bb and cgb; eta_21 use tao_bb and cgb; eta_31 use rho_lb and cgl; eta_41 use tao_lb and cgl
%for eta_13 use tao_13= rho_bl and cgb; eta_23 use tao_bl and cgb; eta_32 use tao_lb and cgl; eta_42 use rho_lb and cgl
%for eta_14 use tao_14= tao_bl and cgb; eta_24 use rho_bl and cgb; eta_34 use tao_ll and cgl; eta_43 use tao_ll and cgl

eta_12= (tao_bb.*cgb)./((omega(i).*L_h).* (2-tao_bb));
eta_13= (rho_bl_lb.*cgb)./((omega(i).*L_h).* (2-rho_bl_lb));
eta_14= (tao_bl_lb.*cgb)./((omega(i).*L_h).* (2-tao_bl_lb));

eta_21= (tao_bb.*cgb)./((omega(i).*L_v).* (2-tao_bb));
eta_23= (tao_bl_lb.*cgb)./((omega(i).*L_v).* (2-tao_bl_lb));
eta_24= (rho_bl_lb.*cgb)./((omega(i).*L_v).* (2-rho_bl_lb));

eta_31= (rho_bl_lb.*cgl)./((omega(i).*L_h).* (2-rho_bl_lb));
eta_32= (tao_bl_lb.*cgl)./((omega(i).*L_h).* (2-tao_bl_lb));
eta_34= (tao_ll.*cgl)./((omega(i).*L_h).* (2-tao_ll));

eta_41= (tao_bl_lb.*cgl)./((omega(i).*L_v).* (2-tao_bl_lb));
eta_42= (rho_bl_lb.*cgl)./((omega(i).*L_v).* (2-rho_bl_lb));
eta_43= (tao_ll.*cgl)./((omega(i).*L_v).* (2-tao_ll));


% Set up the SEA model that includes four subsystem (the updated version on the notebook):

% general formula;
% % P_12= omega(i)* eta_2* M_2* (V_2)^2 + P_21
% % P_d1=omega(i)*eta_1*M_1*(V_1)^2
% 
% %  subsystem 1: bending waves in horizontal beam B_1
% pB_d1= omega(i)*eta_1*M_1*(V_1)^2;
% %  subsystem 2: bending waves in vertical beam B_2
% PB_12= omega(i)* eta_2* M_2* (V_2)^2 + pB_21;
% 
% %subsystem 3: longitudinal waves in horizontal beam B_1
% pL_d1= omega(i)*xi_1*M_1*(V_3)^2;
% %subsystem 4: longitudinal waves in vertical beam B_2
% PL_12= omega(i)* xi_2* M_2* (V_4)^2 + pL_21;


% %% 
% for i=1:length(f)
% om = omega(i)(i);  
% 
% A = [-om*(eta+M1*(eta_12+eta_13+eta_14)), om*M2*eta_21, om*M1*eta_31, om*M2*eta_41;
%     om*M1*eta_13, om*M2*eta_23, -om*(eta+M1*(eta_31+eta_32+eta_34)), om*M2*eta_43;
%     om*M1*eta_12,-om*(eta+M2*(eta_21+eta_23+eta_24)) ,om*M1*eta_32, om*M2*eta_42;
%     om*M1*eta_14, om*M2*eta_24,om*M1*eta_34 ,-om*(eta+M2*(eta_41+eta_42+eta_43))];
% Y = [1;0;0;0];
% v = linsolve(A,Y);
% end



 A = [ -omega(i)*M1*(eta+(eta_12+eta_13+eta_14)),    omega(i)*(M2*eta_21),                       omega(i)*(M1*eta_31),                       omega(i)*(M2*eta_41);
        omega(i)*(M1*eta_13),                        omega(i)*(M2*eta_23),                      -omega(i)*M1*(eta+(eta_31+eta_32+eta_34)),   omega(i)*(M2*eta_43);
        omega(i)*(M1*eta_12),                       -omega(i)*M2*(eta+(eta_21+eta_23+eta_24)) ,  omega(i)*(M1*eta_32),                       omega(i)*(M2*eta_42);  
        omega(i)*(M1*eta_14),                        omega(i)*(M2*eta_24),                       omega(i)*(M1*eta_34) ,                     -omega(i)*M2*(eta+(eta_41+eta_42+eta_43))];
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
rho_ll = 2./(2+6*beta+9*beta.^2); % longitidunal to longitidunal
%%
f1 = figure;
semilogx(f, Lv); 
legend('L_{v1}','L_{v2}','L_{v3}','L_{v4}')
xlabel('f Hz')
ylabel('Velocity Level dB ref. v_1')
grid on;
%%


f2 = figure(Name='Bendıng Coefficients');
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
ylim([0 1])
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
text(beta(index+5).^2,taublRange(fix(length(taublRange)/2)),'\tau_{bb}')

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
f3 = figure;
semilogx(f, v.^0.5); 

legend('v_1','v_2','v_3','v_4')
xlabel('f Hz')
grid on;
%xlim([beta(1)^2,beta(end)^2])
%ylim([0 1])