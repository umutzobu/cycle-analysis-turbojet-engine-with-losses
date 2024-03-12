%!!!For the code to work properly, the file with the "SIunittabless.mat"
%extension must be imported.!!!

clc
clear
fprintf('TURBOJET WITH LOSSES ENGINE CYCLE');
fprintf('\n')
M0 = [0.8, 1, 1.2];
pi_c = 9 ;
T_t4 = 1750 ; %K
hPR = 42800; %kJ/kg
P0 = 29.92;  %kPa
T0 = 229 ; %K
gamma_c = 1.4; 
gamma_t = 1.3;
Cp_c = 1.004 ; %kJ/kgK
Cp_t = 1.239 ; %kJ/kgK
ec = 0.85 ;
et = 0.89;
pi_b = 0.905; 
eta_b = 0.99;
eta_m = 0.98; 
pi_n = 0.98;
P0_P9 = 0.8; 
pi_dmax = 0.95;
ht4=[0 0 0];
f=[0 0.0169 0.0338 0.0507 0.0676]; 
tau_d=1;
load('SIunittabless');
ft=input('initial iterasion guess: ');

%Summary of Equations~Turbojet Engine, CH. 7

for i = 1:length(M0) 
tau_r(i) = 1 + (gamma_c - 1) / 2 .* M0(i).^2;
 if M0(i)<=1
  eta_r(i)=1; 
 else M0(i)>1; 
    eta_r(i)=1-0.075.*(M0(i)-1).^1.35;
end
pi_r(i)=tau_r(i).^(gamma_t/(gamma_t-1));

T = SIunittabless (:, 1);
h3 = SIunittabless (:, 2);

%Diffuser/ Inlet
pi_d(i)=pi_dmax.*eta_r(i); 

%Compressor
tau_c(i)=pi_c.^((gamma_c-1)./(ec*gamma_c));
eta_c(i)=(pi_c.^((gamma_c-1)/(gamma_c)))./(tau_c(i)-1);
T_t3(i)=T0.*tau_r(i).*tau_c(i).*tau_d;
for satir=1:209
    if T_t3(i)==T(satir)
        
    break
    end
end
ht3(i)=h3(satir);
P_t3(i)=P0.*pi_r(i).*pi_c.*pi_d(i);

%Burner
    for b=0:1000 %max iteration number=1000
fi=(ht4(i)-ht3)./((eta_b.*hPR)-ht4(i));
if ft-fi <= 1*10^-4 %Error is 10^-4
    break
end
for j = 1:length(f)-1
    if fi > f(j) && fi < f(j+1)
        break
    end
end
h1 = SIunittabless(167, j); %Our Tt4 temperature is 1750K
h2 = SIunittabless(167, j+1);
ht4(i) = h2-((h2-h1)*(f(j+1)-f(j)))/(f(j+1)-fi);
ft=fi;
    end
    
%Turbine
tau_lambda=(Cp_t*T_t4)/(Cp_c*T0);
tau_t(i)=1-(1./(eta_m.*(1+fi(i))).*(tau_r(i)./tau_lambda).*(tau_c(i)-1));
pi_t(i)=tau_t(i).^(gamma_t*((gamma_t-1)*et));
T_t5(i)=T_t4*(tau_t(i));
P_t5(i)=P0*pi_r(i)*pi_d(i)*pi_c*pi_t(i)*pi_b;

%Nozzle Exit
Pt9_P9(i)=P0_P9*pi_r(i)*pi_d(i)*pi_c*pi_b*pi_t(i)*pi_n;
M9(i)=sqrt((2./(gamma_t-1)).*(Pt9_P9(i).^((gamma_t-1)./gamma_t)));
end

% Draw graphs
figure;

subplot(2, 2, 1);
plot(M0, pi_c, '-o', 'LineWidth', 2);
xlabel('M0');
ylabel('pi_c');
title('Compressor Pressure Ratio');

subplot(2, 2, 2);
plot(M0, pi_t, '-o', 'LineWidth', 2);
xlabel('M0');
ylabel('pi_t');
title('Turbine Pressure Ratio');
subplot(2, 2, 3);
plot(M0, T_t5, '-o', 'LineWidth', 2);
xlabel('M0');
ylabel('T_t5');
title('Turbine Exit Temperature');

subplot(2, 2, 4);
plot(M0, M9, '-o', 'LineWidth', 2);
xlabel('M0');
ylabel('M9');
title('Nozzle Exit Mach Number');

sgtitle('Performance Parameters vs. M0');
