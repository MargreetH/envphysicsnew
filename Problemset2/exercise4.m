%% Load the data
clc; clear all; close all;
importU80V80;


%% a)

U = sqrt(u_80.^2+v_80.^2);
meanU = mean(U);
varU = var(U);

[k,labda,~] = findWeibull(varU,meanU,0.00000001);

u = 0.1:0.01:25;
PDFWeibull = weibullVector(u,labda,k);
plot(u,PDFWeibull);