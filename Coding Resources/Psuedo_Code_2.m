why didnt %% MAE 159 Psuedo Code 2023
% Update 4/24/23
clc;clear all;clf; close all

%% Data Loader %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This would be if you use the method of reading plots from tabular data
% instead of curve fit equations. You data should be stored as a .dat,
% .csv, or .txt file in the working folder.
DeltaMDiv_Cl_Supercritical = table2array(readtable('DeltaMDiv_Cl_Supercritical.csv'))
M_82_Thru_SFC_JT9D = table2array(readtable('82_Thru_SFC_JT9D.csv'))

%% Design Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are parameters you are likely to change less often.
% Your code may or may not use all of the following parameters. You may
% also use many more.
Pax = 275;
W_Cargo = 12000;
R = 4000;         % Range Nautical Miles
TOFL = 6000;
V_Land = 135;     % Landing Approach Speed knots
M = 0.82;          % crusie
Altitude = 35000;
IN = 1.10; % International sizing Choose 1 to turn off
N_Aisle = 2;
N_Abreast = 6;
Fuel_Used = 0.75
N_Engine = 3;
Engine_Fuse = 1
Engine_Wing = 2

%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful atmospehric constants for use in aerodynamic equaitons. Double
% check these numbers using an online atmospheric table for full
% precision.
SpeedSound = 973;
K_Vis = 3.49E-4;  %kinematic viscosity 30,000
d_0 = 9.53E-01;   %density ratio sea level (84F hot day)
d_20000 = 0.5; %20,0000
d_10000 = .7; %10,0000
d_1000 = .9;    %1,000
p = 2.36E-01;     %Pressure ratio (35,0000 ft)

%% Sweep and AR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are assumed constant for now but will need to be later
% explored for some range, ie AR = [6:.5:12].
% This is best done using a 'for' loop arond your final code. For now, keep
% these parameters aside.
AR = 9;
Lam = 35;
%% Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep things neat, labled, and use 'Smart Indent'.

%% Cl Sizing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This loop iterates a cruise Cl for a given configiration.

%% Range Sizing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Aj_f=0; % Ajustment Factor to increase or decrease fuel carried. Starts initially at zero.
Range_Allout=1; % Select a "guess" All Out Range such that the below condition is triggered
R_ao = R + 200 + M*576.48*.75; % All Out Range can be defined here for convience

while abs(R_ao-Range_Allout)>10
    Cl=.5; % Guess Cl
    Cl_f=.1; % Guess Final Cl. This is used to trigger the below conditon. Make sure this is true: abs(Cl_f-Cl)>.005
    while abs(Cl_f-Cl) > .05; % Choose tolerances carefully: Balance speed and accuracy. 'Abs' is nessecary to handle low and high guesses.
        
        % Conditional Statement to Select Supercritical or Conventional airfoil. Omit this for first cut and just use conventional.
        if SupCrit==1
            DeltaM_DIV = interp1(DeltaMDiv_Cl_Supercritical(:,2),DeltaMDiv_Cl_Supercritical(:,1),Cl);
            % Example of interpolaiton from tabulated data (Supercritical).
        else
            DeltaM_DIV= -A*Cl.^2 - B.*Cl + C;
            % Example of using curve fit equation (Conventional).
        end
        M_DIV = (M+.004) - DeltaM_DIV;
        % Conditional Statement to interpolate inbetween sweep angles. Omit this for first cut and just use one angle.
        if SupCrit==1 %Select Supercritical
            if  Lam < 10
                TC = A*M_DIV^2 - B*M_DIV - C;
            elseif(Lam>=10 && Lam<15)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-15)/(20-15)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=15 && Lam<20)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-15)/(20-15)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=20 && Lam<25)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-20)/(25-20)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=25 && Lam<30)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-25)/(30-25)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=30 && Lam<35)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-30)/(35-30)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=35 && Lam<=40)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-35)/(40-35)) + (D*M_DIV^2 - E*M_DIV + F);
            end
        end
        if SupCrit~=1 % Select Conventional airfoil
            if  Lam < 10
                TC = A*M_DIV^2 - B*M_DIV - C;
            elseif(Lam>=10 && Lam<15)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-15)/(20-15)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=15 && Lam<20)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-15)/(20-15)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=20 && Lam<25)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-20)/(25-20)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=25 && Lam<30)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-25)/(30-25)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=30 && Lam<35)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-30)/(35-30)) + (D*M_DIV^2 - E*M_DIV + F);
            elseif(Lam>=35 && Lam<=40)
                TC = ((A*M_DIV^2 - B*M_DIV - C) - (D*M_DIV^2 - E*M_DIV + F)) * ((Lam-35)/(40-35)) + (D*M_DIV^2 - E*M_DIV + F);
            end
        end
        cc=cosd(Lam).^2.*TC.^2.*AR;
        Cl_Takeoff = A.*cc.^2 + B.*cc + C;
        Cl_Landing = A*cc.^3 - B*cc.^2 + C.*cc + D;
        WS_Landing = (V_Land/1.3)^2*((d_0*Cl_Landing)/296);
        R_ao = R + 200 + M*576.48*.75;
        % Conditional Statement to choose weight of fuel for engines JT8D or JT9D. Omit this for first cut and just use JT9D.
        if JT8D==1
            WF_WT = (A*R_ao^2 + B*R_ao+C) + Aj_f
        else
            WF_WT = (A*R_ao^2 + B*R_ao+C) * 0.7820512 + Aj_f % Scaling Factor to change Fuel weight to more efficent JT9D
        end
        WS_Takeoff = WS_Landing/(1-WF_WT*Fuel_Used);
        WS_Cruise = WS_Takeoff*.965;
        Cl_f = WS_Cruise/(1481*M^2*p);
        
        % Conditional Statement to determine whether guess is high or low
        if Cl_f > Cl
            Cl = Cl + 0.01;
        else Cl = Cl - 0.01;
        end
    end
    %% Max Engine Thrust Sizing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Conditional Statement to choose number of Engines.
    % Each equation is the curve for the repsective engine count rearranged to
    % solve for 70% the Max Weight/Thrust ratio as a function of takeoff
    % length.
    if N_Engine ==2
        WT_VT_70 = (d_0 * Cl_Takeoff) * (A * TOFL^2 + B * TOFL - C) / (WS_Takeoff);
    elseif N_Engine ==3
        WT_VT_70 = (d_0 * Cl_Takeoff) * (D * TOFL^2 + E * TOFL - F) / (WS_Takeoff);
    else
        WT_VT_70 = (d_0 * Cl_Takeoff) * (G * TOFL^2 + H * TOFL - I) / (WS_Takeoff);
    end
    
    V_Takeoff = 1.2*((296*WS_Takeoff)/(d_0*Cl_Takeoff))^0.5;
    M_VTakeoff7 = 0.7 * V_Takeoff / 661 / d_0^0.5;
    
    % Conditional Statement to select JT8D or JT9D engines.
    if JT8D==1
        T_0 = 14xxxx; % Max static sea level thrust (Pg. 53)
        % Max Thrust Curve at sea level as function of Mach Takeoff (Pg. 54)
        T_M7 = A*M_VTakeoff7^2 - B*M_VTakeoff7 - C;
        
    else
        T_0 = 45xxxx; % Max static sea level thrust (Pg. 57)
        % Max Dry Thrust Curve at sea level as function of Mach Takeoff (Pg. 61)
        T_M7 = D*M_VTakeoff7^2 - F*M_VTakeoff7 + G;
    end
    
    WT = WT_VT_70*(T_M7/T_0) % Solve for static W/T ratio by scaling W/T @ .70 Velocity
    % by (70% Engine Thrust/Max Engine Thrust)
    %% Weight %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If considering fuselage + wing engines, you need some way of averaging
    % the two engine placement coefficents
    
    % Following Coefficents multipled by some power of the TOGW
    W_Wing = ...
        Length=(3.76*(P/N_Abreast)+33.2) * IN; % Last term international sizing to increase aircraft cabin size
    Diameter=(1.75 * N_Abreast + 1.58 * N_Aisle + 1)*IN;
    W_Fusleage = (0.6727...*Length^0.6*Diameter^ETA*3.75^0.3);
        W_Landing = .04;
    W_Pylon = ...
        W_Total = ...
        W_Power = 1/(WT*3.58);
    W_Fuel = ...
        W_Payload = P*(215)+W_Cargo;
    W_FixE = (132*P+300*N_Engine+260*2+170*(P/50));
    % FC = # Flight Crew
    % CA = # Cabin Atendent
    
    a = W_Total;
    B = W_Fusleage;
    C = W_Landing+W_Pylon+W_Power+W_Fuel+.035-1; % 0.35 Due to one coefficent * TOGW^1
    DD = W_Payload+W_FixE;
    
    W_TO = 100,000; % Guess Weight
    while abs(a*W_TO^1.195 + B*W_TO^.235 + C*W_TO + DD) > 10000 % Ajust tolerances as nessesary; these are obviouly too much
        a*W_TO^1.195 + B*W_TO^.235 + C*W_TO + DD;
        if a*W_TO^1.195 + B*W_TO^.235 + C*W_TO + DD < 10000
            W_TO = W_TO - 1000; % Adjust increment as nessesary to prevent over/undershooting correct value
        else
            W_TO = W_TO + 1000;
        end
    end
    
    % If using composites of advanced tech, reduce the TOGW by each compoment's
    % reduced weight, ie: TOGW - W_Fusleage*TOGW*0.10
    
    S = W_TO/WS_Takeoff;
    b = (AR*S)^5;
    Chord_average = S/b;
    T = W_TO/WT;
    T_Engine = T/N_Engine;
    %% Drag %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Rn = % Reynolds number should be normalized by length, ie.(Re/ft),and taken at 30,000 ft at M = 0.5.
    
    Rn_Wing =Rn*Chord_average;
    Cf_Wing = % Skin friction coefficent of wing
    SWET_Wing =2*1.02*(S-Diameter*30); % Reduce wetted area by that obscured by the fuselage
    K_W = % You may either digitize the plots in 11.3 Shevell or use the formula given in the table in the same book.
    f_Wing = K_W.*Cf_Wing.*SWET_Wing;
    f_Tail = f_Wing*.38;
    
    Rn_Fuselage = Rn*Length;
    Cf_Fuselage = % Skin friction coefficent of wing
    SWET_Fuselage = .9 * pi * Diameter * Length;
    K_f = % Fuselage form factor plot in Shevell
    f_Fuselage = Cf_Fuselage*SWET_Fuselage*K_f;
    
    SWET_Nacel = 2.1*N_Engine*(T_Engine)^.5;
    f_Nacel = 1.25*Cf_Wing*SWET_Nacel;
    f_Pylon = .20*f_Nacel;
    
    f_Total = (f_Wing+f_Fuselage+f_Tail+f_Nacel+f_Pylon) * 1.06;
    C_D0 = f_Total / S;
    e = 1 / (1.035+.38*C_D0*pi*AR);
    %% CLimb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    V_Climb = % Weight_Climb = TOGW * 0.9825 for fuel burned during climb.
    % Atmosphere taken at '20/35 * Cruise altitude'
    % Also TYPO in notes, sample calc V_Cl ~ 485.91
    
    M_Climb = V_Climb * ... % Convert kts to Mach number at '20/35 * Cruise altitude'
        T_Climb = (d_20000 * f_Total * V_Climb ^2)...
        
    % Previous step calculates how much thrust your design aircraft needs.
    
    % Conditional statement to get Max Continuous Climb Thrust and SFC given
    % M_Climb at 20,000 ft
    if JT8D == 1
        T_aJT_D = A*M_Climb^2 - B*M_Climb + c;
        SFC = A*M_Climb^2 + B;
    else
        % Since 20,000 does not exist for JT9D, average between 25k and 15k:
        T_aJT_D = ((A*M_Climb^3 + B*M_Climb^2 - C*M_Climb + D) + (E*M_Climb^2 - F*M_Climb + G)) / 2;
        SFC = ((A*M_Climb^3 + B*M_Climb^2 - C*M_Climb + D) + (E*M_Climb^2 - F*M_Climb + G)) / 2;
    end
    
    T_a = (T_Engine/T_0) * T_aJT_D; % Scale previous values to what design aircraft needs.
    
    Rate_Climb = 101 * V_Climb * (T_a * N_Engine-T_Climb) / (W_TO * .9825);
    
    Time_Climb = Altitude / Rate_Climb;
    Range_Climb = V_Climb * Time_Climb / 60;
    W_FuelClimb = SFC * N_Engine * T_a * Time_Climb / 60;
    
    %% Range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W_0 = W_TO-W_FuelClimb; % Weight at beginning of cruise
    W_1 = (1-WF_WT)* W_TO;
    Cl_Avg = (W_0+W_1)/(...
        C_DI = (Cl_Avg^2)/(...
        C_D = C_D0 + ...
        LD = Cl_Avg/C_D;
    T_Required = (W_0 + W_1...
        T_ReqJT_D =(T_Required/N_Engine)*(T_0/T_Engine);
    
    if JT8D == 1
        SFC_35000 = A*T_ReqJT_D^2 - B*T_ReqJT_D + C;
    else
        SFC_35000 = A*T_ReqJT_D.^3 + B*T_ReqJT_D.^2 - C*T_ReqJT_D;
    end
    
    R_Cruise =(((M.*973.14/1.68781))./SFC_35000).*LD.*log(W_0./W_1); % Note the natural log
    Range_Allout = R_Cruise+Range_Climb;
    WF_WTO = (W_0-W_1)/W_TO;
    
    if Range_Allout < R_ao
        Aj_f = Aj_f + 1; % Ajust tolerances as nessesary; these are obviouly too much
    else
        Aj_f = Aj_f - 1; % Ajust tolerances as nessesary; these are obviouly too much
    end
end
%% Thrust Check and Climb Gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cl_IC = (W_0/S)/(1481*p*M^2);
CDi_IC = (Cl_IC)^2/...
    CD_IC = C_D0...
    LD_IC = Cl_IC/...;
    TReq_IC = W_0/LD_IC/N_Engine;

% Now scale down the above value to that of JTXD:
TReqJT9D_IC = TReq_IC*(T_0/T_Engine);

if JT8D==1
    T_Aval_35000 = % Use the Max continuous cruise line at 35,000 ft to find max cruise thrust
else
    T_Aval_35000 = % Use the Max continuous cruise line at 35,000 ft to find max cruise thrust
end

if TReqJT9D_IC > T_Aval_35000
    fprintf('NOT ENOUGH THRUST TOP OF CLIMB')
end

% NOTE: If you find yourself with insufficent thrust for any of the cruise or climb cases, you will need to
%  decrease the W/T ratio in the 'Max Engine Thrust Sizing' portion. This will
%  be done in a loop similar to that of the Range while loop.

% Climb Gradient Stage 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cmb1C_TO      = Cl_Takeoff/(1.2)^2;
DeltaCmb1_CD0 = ... % Digitize 'Incrimental Profile Drag for High Lift Systems'
Cmb1C_D      = ...
Cmb1LD_TO    = Cmb1C_TO/Cmb1C_D;
Cmb1T_R      = W_TO/Cmb1LD_TO;
Cmb1M        = 1.2.*sqrt(W_TO... % Sea Level Hot Day

if JT8D==1
    Cmb1T_a       = % Use Dry takeoff Engine table, sea level
else
    Cmb1T_a       = % Use Dry takeoff Engine table, sea level
end

GCmb1Grad     =100*((N_Engine-1)...
if N_Engine==2
    if GCmb1Grad<0
        fprintf('C 1 Fail')
    end
elseif N_Engine==3
    if GCmb1Grad<.3
        fprintf('C 1 Fail')
    end
elseif N_Engine==3
    if GCmb1Grad<.5
        fprintf('C 1 Fail')
    end
end

% Stage 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cmb2C_D  = C_D0+DeltaCmb1_CD0...
Cmb2LD_TO = Cmb1C_TO/Cmb2C_D;
Cmb2T_R  = W_TO/Cmb2LD_TO;
GCmb2Grad = 100*((N_Engine-1)...

if N_Engine==2
    if GCmb2Grad<2.4
        fprintf('C 2 Fail')
    end
elseif N_Engine==3
    if GCmb2Grad<2.7
        fprintf('C 2 Fail')
    end
elseif N_Engine==4
    if GCmb2Grad<3
        fprintf('C 2 Fail')
    end
end

% Stage 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cmb3Cl_Max = % Digitize Clean Wing Cl_Max 
Cmb3V      = 1.2*((296*WS_Takeoff)/... Altiude properties at 1000 ft.
Cmb3M      = Cmb3V/...
Cmb3Cl     = Cmb3Cl_Max/(1.2)^2;
Cmb3C_D    = C_D0 + ...
Cmb3LD     = Cmb3Cl/Cmb3C_D;
Cmb3T_R    = W_TO/Cmb3LD ;

if JT8D==1
    Cmb3T_a   = ...
else
    Cmb3T_a   = % Use Max Continuous Max Climb; takeoff, sea level
end

GCmb3Grad  = 100*((N_Engine-1)....
if N_Engine==2
    if GCmb3Grad<1.2
        fprintf('C 3 Fail')
    end
elseif N_Engine==3
    if GCmb3Grad<1.5
        fprintf('C 3 Fail')
    end
elseif N_Engine==4
    if GCmb3Grad<1.7
        fprintf('C 3 Fail')
    end
end

% Approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ApCl       = Cl_Takeoff/1.3^2;
Ap_ClClMax = ApCl/Cl_Takeoff;
ApDeltaCD0 = ... % Digitize 'Incrimental Profile Drag for High Lift Systems'
ApC_D      = C_D0 + ...
ApLD       = ApCl/ApC_D;
ApT_R      = (WS_Landing*S)/ApLD;
Ap_V       = ((296*WS_Landing)/...% Sea Level Hot Day
ApM        = Ap_V /...% Sea Level Hot Day

if JT8D==1
    ApTa      = % Use Dry takeoff Engine table, sea level
else
    ApTa      = % Use Max Continuous % Max Climb; takeoff, sea level
end

GApGrad   =100*(((N_Engine-1)...
if N_Engine==2
    if GApGrad<2.1
        fprintf('Ap Fail')
    end
elseif N_Engine==3
    if GApGrad<2.4
        fprintf('Ap 2 Fail')
    end
elseif N_Engine==4
    if GApGrad<2.7
        fprintf('Ap 2 Fail')
    end
end

%Landing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LCl       = Cl_Landing/1.3^2;
LClClM    = LCl/Cl_Landing ;
LDeltaCD0 = ... % Digitize 'Incrimental Profile Drag for High Lift Systems'
LCD       = C_D0 + ...
LLD       = LCl/LCD;
LT_R      = (WS_Landing*S)/LLD;
LV        = ((296*WS_Landing)/...% Sea Level Hot Day
LM        = LV/...% Sea Level Hot Day

if JT8D==1
    LTa      = % Use Dry takeoff Engine table, sea level
else
    LTa      = % Use Dry takeoff Engine table, sea level
end

GLGrad    =100*(N_Engine*LTa-LT_R)....
if N_Engine==2
    if GLGrad<3.2
        fprintf('Landing Fail')
    end
elseif N_Engine==3
    if GLGrad<3.2
        fprintf('Landing Fail')
    end
elseif N_Engine==4
    if GLGrad<3.2
        fprintf('Landing Fail')
    end
end
%% DOC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = ...
T_CR = ...
T_GM = ...
V_Block = ...
Time_Block = ...
F_Block = ...
% Flight Cost
Passanger = ((165*P+50*P)+W_Cargo)/2000;
Dollar_Hr = ...
CTM_Cr = ...
CTM_Fuel = ...
% Hull Insurace, ect.
Wa= ...

Total_DOC = CTM_Cr + CTM_Fuel + CTM_Hull + CTM_Total + CTM_Depreciation ;
DOC_Passanger = Total_DOC*Passanger/P;
