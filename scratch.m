clc;close all;clear
format long g
pi = 3.141592654;

%% Constants
% Only change this section for different configuration

% C_L loop constants
Conventional_airfoil = 1;
Velocity_approach = 140; % knot
Mach_cruise = 0.82;
Range = 6000; % Nmi
Swept_angle = 35;
AR = 8;

% TOFL constants
Number_of_engine = 3;
Takeoff_field_length = 9000;

% Weight constants
K_w = 1.01; % 1.03 for fuselage engine
Eta = 1.5*2.5; % Ultimate load factor
Constant_Weight_fuselage = 1.1; % 1.1 if 3-class international 1 if not
Taper_ratio = 0.35;
K_f = 11.5; % constant for PAX > 135
PAX = 275; % slang for passenger
N_seats_abreast = 8; % standard
N_aisles = 2; % standard
K_ts = 0.17; % 0.17 for wing engine, 0.25 for fuselage engine
Weight_cargo = 12000; % lb
N_flight_crew = 2;
N_cabin_attendants = 6;

% Drag calculation constants
Fuselage_length = 179; % from tip to tail [ft]
Fuselage_diameter = 20;

% Climb cosntants
Initial_cruise_altitude = 35000; % [ft]

% Thrust check at top climb constant
JT8D = 0; % 1 for JT8D engine and 0 for JT9D engine

% Initalize loop variable
%% C_L loop
% Initalize loop variable
C_L = 0.58;
C_L_final = 0.1;
% While loop boundary and iteration taken from Psuedo_Code.m
while abs(C_L_final-C_L) > .005

    Delta_Mach_div = Delta_Mach_div_function(C_L, Conventional_airfoil);

    % 3. Calculate Divergent Mach number
    Mach_div = (Mach_cruise + 0.004) - Delta_Mach_div;

    % 4. Use Figure 1a to find t/c
    t_c = -0.468*Mach_div + 0.486;

    % 5. Constant: cos^2 t/c AR. Use constant and Fig 3 to find CL_max
    temp = cosd(Swept_angle)^2 * t_c^2 * AR;
    C_L_max_landing = 2.19 + 11.1*temp + -23.2*temp^2;
    C_L_max_takeoff = 1.18 + 12.9*temp + -30.8*temp^2;

    % 6. Calculate wing loading at landing
    sigma = 0.953; % some kind of ratio related to altitude
    WL_landing = (Velocity_approach/1.3)^2*(sigma*C_L_max_landing/296);

    % 7. Crusing velocity and All out range
    Velocity_cruise = round(Mach_cruise * 576.4); % [kts] sqrt(gamma R T) = 576.4
    Range_all_out = Range + 200 + 0.75*Velocity_cruise;

    % 8. Use Figure 4 (Engine JT8D) to find fuel weight to take off weight ratio
    Weight_fuel_takeoff_JT8D = 0.0209 + 1.04E-04*Range_all_out + -5.51E-09*Range_all_out^2;

    % 9. Engine type is JT9D, not JT8D
    SFC_JT9D = 0.61;
    SFC_JT8D = 0.78;
    Weight_fuel_takeoff = Weight_fuel_takeoff_JT8D * (SFC_JT9D/SFC_JT8D) + 0.0257;

    % 10. Take off wing loading
    x = 75/100; % max % fuel at landing
    WL_takeoff =  WL_landing / (1 - x * Weight_fuel_takeoff);

    % 11. Initial crusing wing loading
    WL_initial_crusing = 0.965 * WL_takeoff;

    % 12. Calculate lift coefficient for initial crusing
    C_L_initial_crusing = WL_initial_crusing / (1481 * 0.2360 * Mach_cruise^2);

    C_L_final = C_L_initial_crusing;

    % Conditional Statement to determine whether guess is high or low
    if C_L_final>C_L
        C_L = C_L + 0.01;
    else
        C_L = C_L - 0.01;
    end

end
fprintf('End C_L loop...\n');


%% TOFL

% Weight to Thurst ratio at 0.7 lift off velocity
temp = 31.5*Takeoff_field_length*10^(-3) - 7.45;
Weight_Thrust_0_7_Velocity_liftoff = temp * sigma * C_L_max_takeoff / WL_takeoff;

% 0.7 of Mach number at lift off
Velocity_liftoff = 1.2 * sqrt((296*WL_takeoff) / (sigma*C_L_max_takeoff));
Mach_liftoff = Velocity_liftoff / (661*sqrt(sigma));
Mach_liftoff_0_7 = 0.7 * Mach_liftoff;

% Weight to Thrust ratio
Thrust_JT9D_sea_level_static = 45500; % Sea Level Static Thrust
Thrust_at_0_7_Mach_liftoff = -24567*Mach_liftoff_0_7 + 42600;
Weight_Thrust = Weight_Thrust_0_7_Velocity_liftoff * Thrust_at_0_7_Mach_liftoff / Thrust_JT9D_sea_level_static;

%% Weight
% Weight Wing = Weight_Wing * Weight_takeoff^1.195
Weight_wing = 0.00945 * AR^0.8 * (1 + Taper_ratio)^0.25 * K_w * Eta^0.5 / ( (t_c + 0.03)^0.4 * cosd(Swept_angle) * WL_takeoff^0.695 );

% Weight Fuselage = Weight_fuselage * Weight_takeoff^0.235
l = (3.76*PAX / N_seats_abreast + 33.2) * Constant_Weight_fuselage;
d = (1.75 * N_seats_abreast + 1.58 * N_aisles + 1) * Constant_Weight_fuselage;
Weight_fuselage = 0.6727 * K_f * l^0.6 * d^0.72 * Eta^0.3;

% Weight landing gear = 0.04 * Weight_takeoff
Weight_landing_gear = 0.04;

% Weight nacelle + Weight pylon = Weight_nacelle_pylon * Weight_takeoff
Weight_nacelle_pylon = 0.0555 / Weight_Thrust;

% Weight tail surface = 0.1967 * Weight_wing
Weight_tail_surface = (K_ts + 0.08/Number_of_engine);

% Weight tail surface + wing = Weight_tail_surface_wing * Weight_takeoff^1.195
Weight_tail_surface_wing = (Weight_tail_surface + 1) * Weight_wing;

% Weight power plant = Weight_power_plant * Weight_takeoff
Weight_power_plant = 1/(3.58*Weight_Thrust);

% Weight fuel = Weight_fuel * Weight_takeoff
Weight_fuel = 1.0275 * Weight_fuel_takeoff;

% Weight payload = Weight_payload [lb]
Weight_payload =  215*PAX + Weight_cargo;

% Weight fixed equipment = Weight_fixed_equipment + 0.035*Weight_takeoff
Weight_fixed_equipment = 132 * PAX + 300 * Number_of_engine + 260 * N_flight_crew + 170* N_cabin_attendants;

% Construct Weight polynomial. a*x^1.195 + b*x^0.235 + c*x + d = 0
a = Weight_tail_surface_wing;
b = Weight_fuselage;
c = Weight_landing_gear + Weight_nacelle_pylon + Weight_power_plant + Weight_fuel + 0.035 - 1;
d = Weight_payload + Weight_fixed_equipment;


% Initalize loop variable
Weight_takeoff = 469000;
% While loop boundary and iteration taken from Psuedo_Code.m
while abs(a*Weight_takeoff^1.195 + b*Weight_takeoff^0.235 + c*Weight_takeoff + d) > 100

    % Conditional Statement to determine whether guess is high or low
    if (a*Weight_takeoff^1.195 + b*Weight_takeoff^0.235 + c*Weight_takeoff + d) > 0
        Weight_takeoff = Weight_takeoff + 1000;
    else
        Weight_takeoff = Weight_takeoff - 1000;
    end

end
fprintf('End Weight loop...\n\n');

















function Delta_Mach_div = Delta_Mach_div_function(C_L, Conventional_airfoil)
if Conventional_airfoil == 1
    Delta_Mach_div = -0.348*C_L + 0.191;
elseif Conventional_airfoil == 0
    Delta_Mach_div = -0.179 + 1.07*C_L + -1.84*C_L^2 + 0.873*C_L^3;
else
    fprintf("Airfoil must be either conventional or supercritical.")
end
end
