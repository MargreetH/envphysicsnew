%Computes the cosine of the solar zenith angle at a certain latitude using the julian day
%of the year, at Cabauw, NL. Sets the cosine of the angle to 0 if it is
%smaller than 0.
function cosOfAngle = solarZenithAngle(x)

%For Cabauw
latitude = degtorad(51.965344);
longitude = degtorad(4.897937);
tilt = degtorad(23.45); %tilt of the earth;

xday = floor(x);
xhours = (x - xday)*24;

%calculate declination angle

d = tilt * sin (2*pi / 365 * (284 + xday));

N_d = acos(-tan(latitude)*tan(d)) / pi * 24; %day length in hours

%calculate hour angle
omega = longitude - pi + 2*pi*xhours / 24;

%calculate cosine of solar zenith angle.
cosOfAngle = sin(d)*sin(latitude)+cos(d)*cos(latitude)*cos(omega);

if cosOfAngle < 0
   cosOfAngle = 0; 
end

end