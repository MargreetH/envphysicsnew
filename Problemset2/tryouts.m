clear all;

x = 2:0.01:3;
x1 = x - x(1);
x1 = x1*24;
y = zeros(1,length(x));

for i = 1:1:length(x)
    y(i) = solarZenithAngle(x(i));
    
end


plot(x1,y);
b = y > 0;

dayLength = max(b.*x1) - min(b .* x1)