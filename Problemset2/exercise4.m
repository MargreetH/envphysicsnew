%% Load the data
clc; clear all; close all;
importU80V80;
importUgeoVgeo;

%% a)

U = sqrt(u_80.^2+v_80.^2);
meanU = mean(U);
varU = var(U);

[k,labda,~] = findWeibull(varU,meanU,0.00000001);
delta_u = 0.01;
u1 = 0.1:delta_u:30;
PDFWeibull = weibullVector(u1,labda,k);
plot(u1,PDFWeibull);

%% b)
% Compute annual mean wind power

rho = 1.2;
Rblade = 45;
A_T = pi*Rblade^2;
a= 1/3;

cutInIndex = find(u1==4);
cutOutIndex = find(u1==25);
u1Operational = u1(cutInIndex:cutOutIndex);
PDFWeibullOperational = PDFWeibull(cutInIndex:cutOutIndex);

annualMeanWindPower = sum(2.*rho .* u1Operational.^3 .* A_T .* a.*(1-a).^2 .* PDFWeibullOperational .* delta_u);

fractionTimeOperational = trapz(u1Operational,PDFWeibullOperational);
fractionTimeNotOperational = 1- fractionTimeOperational;

%% compute annual mean theoretic wind power.

Ugeo = sqrt(ugeo.^2+vgeo.^2);
meanUgeo = mean(Ugeo);
varUgeo = var(Ugeo);

[kgeo,labdageo,~] = findWeibull(varUgeo,meanUgeo,0.00000001);
delta_u = 0.01;
u1 = 0.1:delta_u:30;
PDFWeibullgeo = weibullVector(u1,labdageo,kgeo);
plot(u1,PDFWeibullgeo);
PDFWeibullOperationalgeo = PDFWeibullgeo(cutInIndex:cutOutIndex);
annualMeanWindPowergeo = sum(2.*rho .* u1Operational.^3 .* A_T .* a.*(1-a).^2 .* PDFWeibullOperationalgeo .* delta_u);

differenceTheoreticAndActual = abs(annualMeanWindPower - annualMeanWindPowergeo)

