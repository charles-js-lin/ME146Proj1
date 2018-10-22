%% Task 3
% Problem statement:
%       Include effect of solar collector efficiency
%       Determine efficiency of 200 deg. azimuth, tilt angle 36 deg., solar
%       time is between 2-3 PM [14.00-15.00 hrs]
%       
%% Initializing constants
lamb = 37.9;
zeta = 200;
eps = 36;
t = linspace(14.00, 15, 101); %Solar time, [hrs (decimal)]
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
I_d = I_dn.*(cosd(X)*cosd(eps)+ sind(eps)*sind(X).*cosd(E-zeta));

%% Collector Efficiency
ho = 7; %[W/(m^2*K)]
hi = 3.1; %[W/(m^2*K)]
dg = 0.007; %[cm]
kg = 1.3; %[W/mK]
d_ins = 0.06; %[m], from 6 cm
k_ins = 0.018; %[W/mK]
tau = 0.89; %[unitless]
ac = 0.85; %[unitless]
cp = 4186; %[J/(kg*C)]
mdot = 0.0267; %[0.0267 kg/s]

T_amb = 12; %[Celsius]
T_in = 16; %[Celsius]
Area = 1; %assume 1 m^2

COND = 1/(1/(ho*Area) + dg/(kg*Area) + 1/(hi*Area)) + 1/(1/(ho*Area) + d_ins/(k_ins*Area)); %constant
U_loss = COND/Area;
FR = (1-exp(-COND/(mdot*cp)))/(COND/(mdot*cp));

eff = FR*(tau*ac - U_loss./I_d*(T_in-T_amb)); %eff is array

%% Trapezoidal Integration with Efficiency

ID_actual = I_d.*eff; %[array of W/m^2]
ID_total = cumtrapz(t*3600,ID_actual); %[W/m^2] *3600 for [hrs] to [s]
E_total = ID_total*Area; %[W], TOTAL ENERGY

ID_total_noEff = cumtrapz(t*3600,I_d);
E_total_noEff = ID_total_noEff*Area;

plot(t, E_total);
xlabel("Solar time [24-hrs]");
ylabel("Total incident energy [J]");
title("Incident Energy from 2 PM - 3 PM (k_i_n_s = 0.018)");

disp("For k_ins = " + k_ins + "...");
disp("   Total ideal energy (w/out eff.) = " + E_total_noEff(length(E_total_noEff)) + " J");
disp("   Total actual energy = " + E_total(length(E_total)) + " J");
disp("   Efficiency = " + E_total(length(E_total))/E_total_noEff(length(E_total_noEff)));

%% RESULTS

% For k_ins = 0.045...
%    Total ideal energy (w/out eff.) = 3279986.097 J
%    Total actual energy = 2410633.0218 J
%    Efficiency = 0.73495
% 
% For k_ins = 0.018...
%    Total ideal energy (w/out eff.) = 3279986.097 J
%    Total actual energy = 2420375.7269 J
%    Efficiency = 0.73792

%% CONCLUSION

%   After decreasing insulation thermal conductivity to 0.018 W/mK, total
%   energy increased by 9742.71 J and efficiency increased by 0.00297.