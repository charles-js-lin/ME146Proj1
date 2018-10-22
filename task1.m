%% Task 1
% Problem statement:
%     surface azimuth, zeta = 200 deg.
%     tilt angle, eps = 36 deg.
%     day = 4/30/2018
%     time = 1 PM = 13.00 hrs
%     location = Berkeley, CA; lamb = 37.9 deg.
% *All sin/cos/tan and angle entries/calcs are in degrees

%% Initializing constants
lamb = 37.9;
zeta = 200;
eps = 36;
t = 13; %Solar time, [hrs (decimal)]
d = 120;
A = 1310; %[W/m^2]
B = 0.18; %Unitless

%% Intermediate Calculations
a = 15*(t-12); %Hour Angle, [deg]
delta = 23.44*sind(360/365.25*(d-80));
X = acosd(sind(lamb)*sind(delta)+cosd(lamb)*cosd(delta)*cosd(a));

tand_E = sind(a)/(sind(lamb)*cosd(a)-cosd(lamb)*tand(delta));
if a>=0
    if tand_E >= 0
        E = 180 + atand(tand_E);
    elseif tand_E < 0
        E = 360 + atand(tand_E);
    end
elseif a<0
    if tand_E >= 0
        E = atand(tand_E);
    elseif tand_E < 0
        E = 180 + atand(tand_E);
    end
end

%% Final Calculations
% Units should be [W/m^2]

I_dn = A*exp(-B/sind(90-X));
I_d = I_dn*(cosd(X)*cosd(eps)+ sind(eps)*sind(X)*cosd(E-zeta));
disp("Direct solar incident heat flux: " +I_d + " W");

%% Results

% Direct solar incident heat flux: 1048.3376 W