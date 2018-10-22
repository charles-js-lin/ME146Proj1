clc
close all

% constants, inputs, collector's characteristics
A = 1310; % W/m2
B = 0.18;
lambda = 37.9; % Latitude of Berkeley, CA
tau_g = 0.89; %collector glazing transmissivity
alpha_c = 0.85; %collector absorptivity
Cp = 4186; %J/kg*K
h_conv_o = 7; %W/m2*K, outside air convective heat transfer coefficient
delta_g = 0.7*10^-2; %cm converted to m; glazing thickness
Kg = 1.3; %W/m*K; glazing thermal conductivity
h_conv_i = 3.1; %W/m2*K, inside convection coefficient (between absorber and glazing)
delta_ins = 6*10^-2; %cm converted to m; insulation thickness
Kins = 0.045; %W/m*K, insulation thermal conductivity
Ac = 3.25*3; %3 solar collectors, each has area of 3.25 m^2 

% Formulas for efficiency calculation 
% Total conductance is Uloss times area or the collector:
conductance = ((1/(h_conv_o*Ac))+(delta_g/(Kg*Ac))+1/(h_conv_i*Ac))^-1+((1/(h_conv_o*Ac))+(delta_ins/Kins*Ac))^-1;
Uloss = conductance/Ac;

% integration inputs
t0 = 11; % Updated from 10 to 11AM due to shadow being introduced
tfinal = 16; %4PM
h = 0.025; % determined as optimal step size for integration from Task 2
num_steps = (tfinal - t0)/h; % number of steps from step size
T0 = 9+273; % 9C converted to K
Tfinal = 16+273;  % 16C converted to K
Ti = 16+273; % Water inlet temperature, converted from C to K
Tchange = (Tfinal - T0)/num_steps; %use the number of steps to determine change in T for each step

tot_energy = []; %initiate array for Total Energy values
element=1; %initiate the first array index as 1
surf_azimuth_values = [];
surf_inc_values = [];
m_values = [];
m_total = [];
dT_avg = [];
day_num = [];

% Optimal design values determined from previous tasks:
m = 0.0317;
surf_inc = 25;

% Calculate the daily total at two-week intervals, starting with d = 1 and
% finishing at the end of the year.
for d = 1:14:366
outputs = []; %initiate array for Incident Direct Radiation values
index=1; %initiate the first array index as 1
% Ambient temperature varies from 9C to 16C during the time interval
% Start with setting Ta to T0, add linear change in temperature with every
% step. 

Ta = T0; % ambient temperature at the beginning of the process is set to T0
Ta_array = []; % Each temperature value is stored in Ta_array. 
% Last value should be same as Tfinal. 
        
 % intermediate formulas
 alpha = @(t) ((15*(t-12))); % hour angle; a function of t (time of day)
 dec = 23.44*sind((360/365.25)*(d-80)); % declination; varies depending on d (day of the year)
 zen_angle = @(t) rad2deg(acos(sind(lambda)*sind(dec) + cosd(lambda)*cosd(dec)*cosd(alpha(t)))); % zenith angle
 sol_azimuth_tan = @(t) sind(alpha(t))/(sind(lambda)*cosd(alpha(t))-cosd(lambda)*tand(dec)); % tangent of solar azimuth angle
 Int = @(t) A*exp((-B/sind((90)-zen_angle(t)))); % Intensity of direct normal radiation

    for t = t0:h:tfinal %from time t0 to tfinal, with step size h 
    % Determine solar azimuth angle given alpha and tan of solar azimuth angle: 
            if (alpha(t) > 0)
                if(sol_azimuth_tan(t) > 0)
                    sol_azimuth = @(t) 180 + rad2deg(atan(sol_azimuth_tan(t)));
                elseif (sol_azimuth_tan(t) < 0)
                    sol_azimuth = @(t) 360 + rad2deg(atan(sol_azimuth_tan(t)));
                end
            elseif (alpha(t) < 0)
                if(sol_azimuth_tan(t) > 0)
                    sol_azimuth = @(t) rad2deg(atan(sol_azimuth_tan(t)));
                elseif (sol_azimuth_tan(t) < 0)
                    sol_azimuth = @(t) 180 + rad2deg(atan(sol_azimuth_tan(t)));
                end
            end
            
    % Multiplied by "clear conditions" percentage (46% of hours are
    % cloudy, so total output should be expressed as  54% of what 
    % was originally calculated).
    surf_azimuth = sol_azimuth;
    Id = @(t) 0.54*Int(t)*(cosd(zen_angle(t))*cosd(surf_inc)+sind(surf_inc)*sind(zen_angle(t))*cosd(sol_azimuth(t)-surf_azimuth(t)));
    %surf_azimuth = 214;
    %Id = @(t) 0.54*Int(t)*(cosd(zen_angle(t))*cosd(surf_inc)+sind(surf_inc)*sind(zen_angle(t))*cosd(sol_azimuth(t)-surf_azimuth));
    % Heat removal factor:
    Fr = (1 - exp(-conductance/(m*Cp)))/(conductance/(m*Cp));
    efficiency =  @(t) Fr*((tau_g*alpha_c) - (Uloss/Id(t))*(Ti - Ta));
    % Multiplied by Ac to get the value for the entire collector area
    outputs(index) = Id(t)*efficiency(t)*Ac;
    Ta_array(index) = Ta;
    index = index + 1;
    Ta = Ta + Tchange;
    end
    % Use trapz MATLAB function, using time_steps as x and outputs(t) as y
    time_steps = [t0:h:tfinal]*3600; % Convert time steps from hours to 
    % seconds for integration purposes (since W = J/s)
    tot_energy (element) = trapz(time_steps,outputs);
    day_num (element) = d;
    element = element + 1;
end

tot_energy_transposed = tot_energy';
day_num_transposed = day_num';
% Concatenate the column vectors for total energy, azimuth, and inclination
% into a results matrix
results = (horzcat(tot_energy_transposed,day_num_transposed));

% Integrate over the entire year to find the amount of energy collected for
% the year
yearly_energy = trapz(day_num,tot_energy);
format long g
disp(yearly_energy); 

% Part (a): Total heat collected
% 19546821131.4158 J = 19.55 GJ with original design
% 20463180561.8825 J = 20.46 GJ with new disign

% Part (b): Percentage increase in yearly energy over origibal design
% Improvement: 916359430.466702 J = 916.36 MJ
% Percent increase: 0.04688 = 4.7% improvement

% Part (c): How extra energy from motorized frame design compares to shadow loss
% Value before shadow was introduced: 23028578254.8 J = 23.03 GJ
% Value after shadow: 19.55 GJ
% Value after shadow with improved design: 20.46 GJ
% Loss from shadow: 3481757123.3736 J = 3.48 GJ
% Gain from new design: 916.36 MJ
% Percentage of shadow loss compensated by new design: 26.3%
% (916359430.466702/3481757123.3736 = 26.3%)

% The motorized frame collector would not produce enough extra
% energy to compensate for the loss due to new building's shadow.
% The new design would only be able to offset a quarter of the energy lost
% due to the shadow being introduced.


