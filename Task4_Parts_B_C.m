clc
clear all
close all

% constants, inputs, collector's characteristics
A = 1310; % W/m2
B = 0.18;
lambda = 37.9; % Latitude of Berkeley, CA
d = 120; % Day of the year; Day 120 is Apr 30, 2018
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

% Mass flow rate: Optimize to maximize volume of water heated.
%m = 0.0614; %mass flow rate, in kg/s
% Heat removal factor:
%Fr = (1 - exp(-conductance/(m*Cp)))/(conductance/(m*Cp));

% intermediate formulas
alpha = @(t) ((15*(t-12))); 
dec = 23.44*sind((360/365.25)*(d-80)); %declination
zen_angle = @(t) rad2deg(acos(sind(lambda)*sind(dec) + cosd(lambda)*cosd(dec)*cosd(alpha(t))));
sol_azimuth_tan = @(t) sind(alpha(t))/(sind(lambda)*cosd(alpha(t))-cosd(lambda)*tand(dec)); % tangent of solar azimuth
Int = @(t) A*exp((-B/sind((90)-zen_angle(t)))); %Intensity

% integration inputs
t0 = 10; % 10AM
tfinal = 16; %4PM
h = 0.025;
num_steps = (tfinal - t0)/h; %number of steps given step size
Ti = 16+273; % Water inlet temperature, converted from C to K
% Ambient temperature varies from 9C to 16C during the time interval
T0 = 9+273; % 9C converted to K
Tfinal = 16+273;  % 16C converted to K
Tchange = (Tfinal - T0)/num_steps; %use the number of steps to determine change in T for each step

tot_energy = []; %initiate array for Total Energy values
element=1; %initiate the first array index as 1
surf_azimuth_values = [];
surf_inc_values = [];
m_values = [];
m_total = [];
dT_avg = [];

m = 0.0317;
%optimal mass flow: 0.0317 kg/sec, achieving 684.72kg of heated water in 
%the tank, for the average temperature increase of 49.12C (avg temperature 
%of water in the tank at the end of the day is 65.12C, assuming the tank
%is perfectly insulated.
%for surf_azimuth = 90:1:270
   surf_azimuth = 214;
   %for surf_inc = 20:1:60
   surf_inc = 25;
        outputs = []; %initiate array for Incident Direct Radiation values
        index=1; %initiate the first array index as 1
        Ta = T0; % ambient temperature at the beginning of the process is set to T0
        Ta_array = [];
        dT = [];
        dT_sum = 0;
        for t = t0:h:tfinal %from time t0 to tfinal, in 0.025 hr increments
        % Rule to determine solar azimuth angle given alpha and tan of solar
        % azimuth angle. 
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
        % Multiplied by Ac to get the value for the entire collector area
            Id = @(t) Ac*Int(t)*(cosd(zen_angle(t))*cosd(surf_inc)+sind(surf_inc)*sind(zen_angle(t))*cosd(sol_azimuth(t)-surf_azimuth));
            % Heat removal factor:
            Fr = (1 - exp(-conductance/(m*Cp)))/(conductance/(m*Cp));
            efficiency =  @(t) Fr*((tau_g*alpha_c) - (Uloss/Id(t))*(Ti - Ta));
            outputs(index) = Id(t)*efficiency(t);
            dT (index) = (outputs(index)/(m*Cp));
            dT_sum = dT_sum + dT(index);
            Ta_array(index) = Ta;
            index = index + 1;
            Ta = Ta + Tchange;
        end
        % Use trapz MATLAB function, using time_steps as x and outputs(t) as y
        time_steps = [t0:h:tfinal]*3600;
        tot_energy (element) = trapz(time_steps,outputs);
        %surf_azimuth_values (element) = surf_azimuth;
        %surf_inc_values (element) = surf_inc;
        m_values (element) = m;
        m_total (element) = m*(tfinal-t0)*3600;
        dT_avg (element) = dT_sum/num_steps;
        element = element + 1;
    %end
%end
%end

% Transpose vectors for total energy, azimuth, and inclination to columns
    tot_energy_transposed = tot_energy';
    %surf_azimuth_transposed = surf_azimuth_values';
    %surf_inc_transposed = surf_inc_values';
    m_values_transposed = m_values';
    m_total_transposed = m_total';
    dT_avg_transposed = dT_avg';
    
% Concatenate the column vectors for total energy, azimuth, and inclination
% into a results matrix
    results = (horzcat(tot_energy_transposed,m_values_transposed,m_total_transposed,dT_avg_transposed));

% Find the maximum value of the first column of the results matrix
% (highest energy produced), and display the corresponding row of the 
% matrix (energy value and the corresponding azimuth and inclination).
    [~,x] = max(results(:,3));
    format long g
    answer = results(x,:);
    %disp(answer);
    %disp(results);

%PART B
Tstart = 16;
Tfinish = 65;
m_req = 350*3.79; % converting 350 gallons of water to kilograms   
m_left = m_req - m_total; % how much water still needs to be heated by CH4
Q = m_left*Cp*(Tfinish-Tstart); % heat required to heat the remaining water
CH4_heating = 50050000; % in J/kg; lower heating value of CH4
CH4_kg_ideal = Q/CH4_heating; % how much CH4 it would take if the process was 100% efficient
CH4_kg_90eff = CH4_kg_ideal/0.9; % how much CH4 it would take for the 90% efficient process

disp("-----PART B-----");
disp("Total amount of water to be heated:   " + m_req + " kg");
disp("Amount of water heated by collector energy:   " + m_total + " kg");
disp("Amount of water that remains to be heated by natural gas:   " + m_left + " kg");
disp("Amount of heat required to heat the remaining water:   " + Q + " J");
disp("Natural gas required per year to heat remaining water:   " + CH4_kg_90eff + " kg");
%answer: 2.9224kg of natural gas at 90% efficiency would be required to heat the remaining water.

%PART C
num_days = 200;
m_yearly = num_days*m_total; % total mass of water heated in 200 days
% Heat required is mass of water times Cp times change in temperature:
Q_yearly = m_yearly*Cp*(Tfinish-Tstart);
% Calculating how much CH4 would be required to generate this heat
% (assuming 90% efficiency):
CH4kg_saved = (Q_yearly/CH4_heating)/0.9;%answer:623.5793kg
% Molar mass of CH4 and CO2 required for these calculations:
M_CH4 = 16; %grams per mol
M_CO2 = 44; %grams per mol
% Number of moles is equal to mass divided by molar mass:
mol_CH4 = (CH4kg_saved*1000)/M_CH4; %moles, with mass converted to grams first
%number of moles is the same, from coefficients of chemical reaction:
% CH4 + O2 = CO2 + 2H2O
mol_CO2 = mol_CH4; 
% Mass of CO2 is equal to molar mass of CO2 times the number of moles:
CO2kg_saved = (M_CO2*mol_CO2)/1000;%answer:1714.8432kg

disp("-----PART C-----");
disp("Amount of CO2 that would be prevented from being released yearly:   " + CO2kg_saved + " kg");
