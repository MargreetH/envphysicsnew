%% refreshes the data
clc; clear all; close all;
importyear; %Simply imports the data cont
importSolar;
%Remove first nonsensical element

%% compute mean for year 2012

currentYear = hyea(1);
counter = 1;
maxIterations = 10^6;
cond = true;
wrongValue = -9999;
i = 0;
correctValues = 0;
currentLWupSum = 0;
currentLWdnSum = 0;
currentSWupSum = 0;
currentSWdnSum = 0;
currentLWdnTopAtmosSum = 0;
currentSWdnTopAtmosSum = 0;

numberOfYears = countYears(hyea);
averages = zeros(4,numberOfYears);
yearCounter = 1;
clc;
while cond
    
  i = i + 1; %Incremented each loop cycle;
  if (hyea(i) ~= currentYear) 
      averages(1,yearCounter) = currentLWupSum /correctValues;
      averages(2,yearCounter) = currentLWdnSum /correctValues;
      averages(3,yearCounter) = currentSWupSum /correctValues;
      averages(4,yearCounter) = currentSWdnSum /correctValues;      
      currentYear = hyea(i);
      correctvalues = 0;
      yearCounter = yearCounter+1;
  end
  
  %Check for incorrect data
  cond1 = LW_up(i) ~= wrongValue;
  cond2 = LW_dn(i) ~= wrongValue;
  cond3 = SW_up(i) ~= wrongValue;
  cond4 = SW_dn(i) ~= wrongValue;
  
  if (cond1 && cond2 && cond3 && cond4) %Datapoint is valid
      correctValues = correctValues + 1;
      currentLWupSum = currentLWupSum + LW_up(i);
      currentLWdnSum = currentLWdnSum + LW_dn(i);
      currentSWupSum = currentSWupSum + SW_up(i);
      currentSWdnSum = currentSWdnSum + SW_dn(i);
  else %Datapoint is invalid
      %do nothing
  end
  
  if numberOfYears == 1 && i == length(hyea) %we have only one year, and 
      averages(1) = currentLWupSum /correctValues;
      averages(2) = currentLWdnSum /correctValues;
      averages(3) = currentSWupSum /correctValues;
      averages(4) = currentSWdnSum /correctValues;
      cond = false;
  end
    
 if i == length(hyea)
     cond = false;
 end

end

reflection = averages(3)/averages(4)


%% annual mean top atmosphere at cabauw

x = 1:0.01:365; %In julian days
y = zeros(1,length(x));
S_0 = 1360; 

for i = 1:1:length(x)
   y(i)=solarZenithAngle(x(i)) * S_0;
end

average = mean(y)
%%