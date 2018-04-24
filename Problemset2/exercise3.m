%% constants
clc;
clear all;
close all;
T_sfc = 293;
T_sfc2 = 283;
p_sfc = 1.018 * 10^5;
R_d = 287.04;
rho_sfc =  p_sfc / (R_d * T_sfc);
g=9.81;


%% a)
dz =1;
z = 0:dz:10^4;

p_isoT = p_sfc .* exp(-g .* z/(R_d * T_sfc));


figure;
plot(z,p_isoT,'LineWidth',2);
xlabel('z [m]')
ylabel('p [Pa]')
title('Pressure decline with height for an isothermal atmosphere')

set(gca,'FontSize',10) % make fontsize bigger
set(gcf,'color','w'); % Set bg color to white





rho_isoT = p_isoT ./(R_d*T_sfc);
figure;
plot(z,rho_isoT,'LineWidth',2);
xlabel('z [m]')
ylabel('density [kg/m^3]')
title('Density decline with height for an isothermal atmosphere')

disp('pressure at 5000 m');
disp(p_isoT(5001));
disp('density at 5000 m');
disp(rho_isoT(5001));

%% b) do the same

p_isoT2 = p_sfc .* exp(-g .* z/(R_d * T_sfc2));
rho_isoT2 = p_isoT2 ./(R_d*T_sfc2);



figure;
plot(z,p_isoT, z, p_isoT2,'LineWidth',2);
xlabel('z [m]')
ylabel('p [Pa]')
legend('Tsfc = 293', 'Tsfc = 283')
title('Pressure decline with height for an isothermal atmosphere')
set(gca,'FontSize',10) % make fontsize bigger
set(gcf,'color','w'); % Set bg color to white

figure;
plot(z,rho_isoT, z,rho_isoT2,'LineWidth',2);
xlabel('z [m]')
ylabel('density [kg/m^3]')
legend('Tsfc = 293', 'Tsfc = 283')
title('Density decline with height for an isothermal atmosphere')

set(gca,'FontSize',10) % make fontsize bigger
set(gcf,'color','w'); % Set bg color to white
horzPressureGradient = (p_isoT2(5001)-p_isoT(5000))/(500*10^3);

rhoInBetween = (rho_isoT(5001) + rho_isoT2(5001))/2;
du_dt = -horzPressureGradient /rhoInBetween

%% Numerical integration

%Create the z vector to hold the steps
dZ = 0.001; %can be changed for finetuning;
totalZ = 10^4;
zvector = 0:dZ:totalZ;
gamma = 6E-3;
Tk = T_sfc;
Tk_plus_one = 0;

%Create vectors to hold the values for the pressure and density
p_Tvariable = zeros(1,length(zvector));
rho_Tvariable = zeros(1,length(zvector));
p_Tvariable(1) = p_sfc;
rho_Tvariable(1) = rho_sfc;
T_var = zeros(1,length(zvector));
T_var(1) = T_sfc;

for i=1:1:length(zvector)-1
    zi = zvector(i);
    T_var(i+1) = 2 * (T_var(1) - gamma * (zi + 0.5*dZ))-T_var(i);

    
    p_Tvariable(i+1) = p_Tvariable(i) .* exp(-g .* dZ/(R_d * T_var(i)));
    rho_Tvariable(i+1) =  p_Tvariable(i+1)/(R_d * T_var(i));
    
end


figure;
plot(z,p_isoT, zvector, p_Tvariable,'LineWidth',2);
xlabel('z [m]')
ylabel('p [Pa]')
legend('Isothermal atmosphere', 'Varying temperature, gamma = 6E-3')
title('Pressure decline with height for an isothermal atmosphere, Tsfc = 293')
set(gca,'FontSize',10) % make fontsize bigger
set(gcf,'color','w'); % Set bg color to white

figure;
plot(z,rho_isoT, zvector,rho_Tvariable,'LineWidth',2);
xlabel('z [m]')
ylabel('density [kg/m^3]')
legend('Isothermal atmosphere', 'Varying temperature, gamma = 6E-3')
title('Density decline with height for an isothermal atmosphere, Tsfc = 293')

set(gca,'FontSize',10) % make fontsize bigger
set(gcf,'color','w'); % Set bg color to white

%Code to calculate difference
if length(zvector) > length(z)
    shortV = z;
    longV = zvector;
    shortP = p_isoT;
    longP = p_Tvariable;
else  
    shortV = zvector;
    longV = z;
    longP = p_isoT;
    shortP = p_Tvariable;
end
difference = 0;

for i=1:1:length(shortV)
   index = find(longV == shortV(i)); 
   difference = abs(longP(index)-shortP(i));
   if difference > 0.05
       disp('difference larger than 5% at z=');
       disp(shortV(i));
       break;
   end   
end

%% d)

%compute q_v
e_sat = 610.78 .* exp(17.2694 .* (T_var-273.16)./(T_var-35.86));
r_sat = 0.622 * e_sat ./ p_Tvariable;
m_v = rho_Tvariable .* r_sat;
q_v = m_v ./ (m_v + rho_Tvariable);

%Plot q_v versus z
figure;
plot(zvector,q_v,'LineWidth',2);
xlabel('z [m]')
ylabel('q_v')
title('Specific humidity as a function of height')
set(gca,'FontSize',10) % make fontsize bigger
set(gcf,'color','w'); % Set bg color to white

%Compute water vapor path
WVP_1 = dZ/2 * (rho_Tvariable(1)*q_v(1) + rho_Tvariable(end)*q_v(end));
for i=2:1:(length(zvector)-1)
    WVP_1 = WVP_1 + dZ * (rho_Tvariable(i)*q_v(i));
end

display('WVP: ');
display(WVP_1);

%% d ii)
dZ = 1; %can be changed for finetuning;
totalZ = 10^4;
zvector = 0:dZ:totalZ;

temperatures = 260:1:310;
gamma = 6E-3;
%Create the z vector to hold the steps

%create new vectors to hold the new solutions for each T
p_Tvariableii = zeros(length(temperatures),length(zvector));
rho_Tvariableii = zeros(length(temperatures),length(zvector));
T_varii = zeros(length(temperatures),length(zvector));
WVP = zeros(1,length(temperatures));

for j=1:1:length(temperatures)
    
    %Create vectors to hold the values for the pressure and density
    p_Tvariableii(j,1) = p_sfc;
    rho_Tvariableii(j,1) = rho_sfc;
    T_varii(j,1) = temperatures(j);

    for i=1:1:(length(zvector)-1)
        zi = zvector(i);
        T_varii(j,i+1) = 2 * (T_varii(j,1) - gamma * (zi + 0.5*dZ))-T_varii(j,i);
        p_Tvariableii(j,i+1) = p_Tvariableii(j,i) .* exp(-g .* dZ/(R_d * T_varii(j,i)));
        rho_Tvariableii(j,i+1) =  p_Tvariableii(j,i+1)/(R_d * T_varii(j,i));
        
    end
    
    %compute q_v
    e_sat = 610.78 .* exp(17.2694 .* (T_varii(j,:)-273.16)./(T_varii(j,:)-35.86));
    r_sat = 0.622 * e_sat ./ p_Tvariableii(j,:);
    m_v = rho_Tvariableii(j,:) .* r_sat;
    q_v2 = m_v ./ (m_v + rho_Tvariableii(j,:));   
    WVP(j) = dZ/2 * (rho_Tvariableii(j,1)*q_v2(1) + rho_Tvariableii(j,end)*q_v2(end));   

for k=2:1:(length(zvector)-1)
    WVP(j) = WVP(j) + dZ * (rho_Tvariableii(j,k)*q_v2(k));
end

end

figure;
plot(temperatures,WVP,'LineWidth',2);
xlabel('surface temperature [K]')
ylabel('Water Vapor Path [kg/m^3]')
title('Maximum Water Vapor Path as a function of surface temperature')
set(gca,'FontSize',10) % make fontsize bigger
set(gcf,'color','w'); % Set bg color to white






















