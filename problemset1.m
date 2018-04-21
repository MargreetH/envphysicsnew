
%% c)
clc;
clear all;
close all;

%All constants
N_A = 6.022*10^23; % per mol

%Number densities in original spent fuel mixture. per cm^3
N_FP_PF = 1.93*10^21;
N_235_PF = 1.45*10^21;
N_238_PF = 4.59*10^22;
N_TOTAL = N_FP_PF + N_235_PF + N_238_PF;
rho_U = 19.1; %g/cm^3
barnConversion = 10^(-24); % to cm^2

%Constants from table
nu235 = 2.43;
 %in barns:
 sigmaF235 = 584 * barnConversion;
 sigmaC235 = 98* barnConversion;
 sigmaC238 = 2.7* barnConversion;
 sigmaCFP = 135* barnConversion;
 sigmaCH2O = 0.68* barnConversion;
 sigmaCHNO3 = 2.22* barnConversion;
 M_H2O = 18; % g per mol
 M_U = 237.9004; % g per mol.
 concentration_HNO3 = 0.007; % mol per cubic centimeter
 
 %Calculate fixed number densities for H2O and NHO3:
 N_H2O = N_A/M_H2O; % number per cm^-3
 N_HNO3 = concentration_HNO3*N_A; % number per cm^-3
 
 %Vector containing the different concentrations for the SF mixture:
 concentrationSF = 0:0.01:20; %g/cm^-3;
 kINF = zeros(1,length(concentrationSF));
 
 epsilon = 1;
 f = 1;
 p = 1;
 
 disp(M_U * N_A);
 
 for i = 1:1:length(concentrationSF)
     
     %Calculate all number densities:
     N_TOTAL_DILUTED = concentrationSF(i) / M_U * N_A;
     N_FP_LM = N_FP_PF/ N_TOTAL * N_TOTAL_DILUTED;
     N_235_LM = N_235_PF /N_TOTAL * N_TOTAL_DILUTED;
     N_238_LM = N_238_PF /N_TOTAL * N_TOTAL_DILUTED;
     
     sigma_TOTAL_FISSION = sigmaF235 * N_235_LM; %Only Uranium 235 contributes
     sigma_TOTAL_ABSORPTION = (sigmaC235+sigmaF235) * N_235_LM + sigmaC238 * N_238_LM + sigmaCFP * N_FP_LM + sigmaCH2O * N_H2O + sigmaCHNO3 * N_HNO3;
     
     nu = nu235 * sigma_TOTAL_FISSION / sigma_TOTAL_ABSORPTION;
     
     kINF(i) = nu * epsilon * f * p;
 end
 
 checkittt = 105.7 .* concentrationSF ./ (32.1 + 70.6 .* concentrationSF );

plot(concentrationSF, kINF,'LineWidth',2);

set(gca,'FontSize',10) % make fontsize bigger
set(gcf,'color','w'); % Set bg color to white

%Polish the axes
newLim = get(gca,'XLim'); 
newx = linspace(newLim(1), newLim(2), 20); 



xlabel('spent fuel density in g/cm^3')
ylabel('infinite multiplication factor k')
title('Infinite multiplication factor as a function of the spent fuel density')

%% PART 2

R = 1.3;
H = 3;
D = 1.3;
alpha0 = 2.4048;

B_g = sqrt( (pi/H)^2 + (alpha0/R)^2) * 100; % 1/cm

rholmsf = 0.80;
%disp(D * B_g^2)

%k_for_e = 105.7 * rholmsf / (70.6 * rholmsf + 32.1 + D * B_g^2 );


%RETRY, hopefully correct:

k_inf_e = 143.5 * rholmsf / (43.61 + 96.81 * rholmsf);
L_squared = D / (0.0706 * rholmsf+0.0321);
P_for_e = (1/ L_squared * B_g^2);

k_for_e = k_inf_e * P_for_e;



%% PART 3 calculating l_eff

R = 1.3*10^(-2);
H = 3*10^(-2);
D= 1.3;
alpha0 = 2.4048;

B_g = sqrt( (pi/H)^2 + (alpha0/R)^2);

rho_SF_mistake = 0.95; %g/cm^3
t_d = 13; %s
u = 2*10^7; %m/s
beta = 0.006;



sigma_a = 0.0706 * rho_SF_mistake + 0.0321;
L = sqrt(D/sigma_a);

l = 1 / ( sigma_a * u * (1 + L^2 * B_g ^ 2));
l_eff = (1-beta) * l + beta * (l + t_d);

k_inf_h = 143.5 * rho_SF_mistake / (43.61 + 96.81 * rho_SF_mistake);
P_for_h = (1/ L_squared * B_g^2);

k_for_h = k_inf_h * P_for_h;


waiting_time = log(10) * l_eff / (k_for_h - 1)
