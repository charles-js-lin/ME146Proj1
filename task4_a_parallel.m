%% Task 4 - PARALLEL
% Problem statement:
%   * Design params: zeta (surf. az); eps (surf. tilt); mdot (mass flow
%           rate)
%   * T_out: 65 oC
%   * Solar time: 10.00 -> 16.00 (10 AM - 4 PM)
%   * qdot [W] = mdot*cp*dT = E_tot
%       * cp*dT = 4186*(65-16) = 205114
%       * mdot = E_tot/205114

%% Initializing constants
% zeta and eps will be initialized before solar flux calculations

lamb = 37.9;
t = linspace(10.00, 16, 601); %Solar time, [hrs (decimal)]
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

%% Incident Direct Solar Flux
% Units should be [W/m^2]

I_dn = A*exp(-B./sind(90-X));
eps = linspace(0,90,601); %DESIGN PARAM: SURFACE TILT ANGLE
zeta = linspace(0,360,601); %DESIGN PARAM: SURFACE AZIMUTH ANGLE

maxI_d = 0;
OPT_eps = 0;
OPT_zeta = 0;

for i=1:601
    for j=1:601
        I_d = I_dn.*(cosd(X).*cosd(eps(i))+ sind(eps(i)).*sind(X).*cosd(E-zeta(j)));
        if max(I_d)>maxI_d
            maxI_d = max(I_d);
            OPT_eps = eps(i);
            OPT_zeta = zeta(j);
        end
    end
end

% Maximized I_d
I_d = I_dn.*(cosd(X).*cosd(OPT_eps)+ sind(OPT_eps).*sind(X).*cosd(E-OPT_zeta));

%% Collector Efficiency
ho = 7; %[W/(m^2*K)]
hi = 3.1; %[W/(m^2*K)]
dg = 0.007; %[cm]
kg = 1.3; %[W/mK]
d_ins = 0.06; %[m], from 6 cm
k_ins = 0.045; %[W/mK]
tau = 0.89; %[unitless]
ac = 0.85; %[unitless]
cp = 4186; %[J/(kg*C)]
mdot = 0.0317/3; %DESIGN PARAM: MASS FLOW RATE, MUST BE DIVIDED BY 3

T_amb = linspace(9, 16, 601); %[Celsius] %VARIES LINEARLY between 9-16 oC
T_in = 16; %[Celsius]
Area = 3.25; %assume 3.25 m^2, ONE COLLECTOR

COND = 1/(1/(ho*Area) + dg/(kg*Area) + 1/(hi*Area)) + 1/(1/(ho*Area) + d_ins/(k_ins*Area)); %constant
U_loss = COND/Area;
FR = (1-exp(-COND/(mdot*cp)))/(COND/(mdot*cp));

eff = FR*(tau*ac - U_loss./I_d.*(T_in-T_amb)); %eff is array, includes I_d

%% Trapezoidal Integration with Efficiency

ID_actual = I_d.*eff; %[array of W/m^2]
ID_total = cumtrapz(t*3600,ID_actual); %[J/m^2], *3600 to convert from [hrs] to [s]
E_coll = ID_total*Area; %[J], TOTAL ENERGY
E_total = 3*E_coll; % PARALLEL - 3*individual collectors

%% Results

plot(t, E_total);
xlabel("Solar time [24-hrs]");
ylabel("Total incident energy [J]");
title("Incident Energy from 10 AM - 4 PM");

disp("Parallel configuration: ");
disp("   Total energy = " + E_total(length(E_total)) + " J");
disp("   Max efficiency = " + max(eff));
disp("   Min efficiency = " + min(eff));
disp("   Optimal surface tilt angle epsilon: " + OPT_eps);
disp("   Optimal surface azimuth angle zeta: " + OPT_zeta);
