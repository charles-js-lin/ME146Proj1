%% Task 2
% Problem statement:
%     surface azimuth, zeta = 200 deg.
%     tilt angle, eps = 36 deg.
%     day = 4/30/2018
%     time =10 AM to 4:00PM = 10.00 - 16.00 hrs
% *All sin/cos/tan and angle entries/calcs are in degrees

%% Initializing constants
lamb = 37.9;
zeta = 200;
eps = 36;
t = linspace(10.00, 16.00, 601); %Solar time, [hrs (decimal)]
d = 120;
A = 1310; %[W/m^2]
B = 0.18; %Unitless

%% Intermediate Calculations
a = 15*(t-12); %Hour Angle, [deg] (array)
delta = 23.44*sind(360/365.25*(d-80));
X = acosd(sind(lamb)*sind(delta)+cosd(lamb)*cosd(delta)*cosd(a));

for i = 1:length(a)
    tand_E(i) = sind(a(i))/(sind(lamb)*cosd(a(i))-cosd(lamb)*tand(delta));
    if a(i)>=0
        if tand_E(i) >= 0
            E(i) = 180 + atand(tand_E(i));
        elseif tand_E(i) < 0
            E(i) = 360 + atand(tand_E(i));
        end
    elseif a(i)<0
        if tand_E(i) >= 0
            E(i) = atand(tand_E(i));
        elseif tand_E(i) <= 0
            E(i) = 180 + atand(tand_E(i));
        end
    end
end

%% Final Calculations
% Units should be [W/m^2]

I_dn = A*exp(-B./sind(90-X));
I_d = I_dn.*(cosd(X)*cosd(eps)+ sind(eps)*sind(X).*cosd(E-zeta));

%% Time Step Verification and Integration
for j = 1:length(I_d)-1
    new(j)=abs(100*(I_d(j+1)-I_d(j))/I_d(j+1));
end

% Making a table of dI_d over dt
for i=1:length(t)-1
    dt_arr(i) = t(i+1)-t(1);
end
dI_d = new';
dt = dt_arr';
table = [dt dI_d];


Z = cumtrapz(t*3600,I_d); %Converting t from [hrs] to [s]
plot(t,Z);
xlabel("Solar time [24-hrs]");
ylabel("Total incident energy [J]");
title("Incident Energy from 10 AM - 4 PM");

disp("Total energy: " + Z(length(Z)) + " J");

%% Results

% Total incident energy: 20036071.2047 J