%%% LPG ACTUATOR MODEL
% Auth: Thomas W. C. Carlson
% clc, clear all
% Parameters given by subsystem
%% Spool Input Subsystem
% No direct parameters in this subsystem
% To design a signal open the Spool Input subsystem and open the Signal
% Editor, then launch the Signal Editor UI and create the .mat scenario
% desired.
% TODO: Script the writing of .mat scenarios based on commanded operations
%% Exhaust R Subsystem
% Reservoir Parameters
ExhRR_Pres = 1;                         % [atm], Reservoir Pressure
ExhRR_Temp = 273.15;                    % [K], Reservoir Temperature
ExhRR_Diam = 0.125;                     % [in], Exhaust Hose Inlet Diameter
ExhRR_Area = (ExhRR_Diam/2)^2*pi;       % [in^2], Exhaust Cross-section Area
% Pipe Parameters
ExhPR_Length = 1;                       % [in], Exhaust Pipe Length
ExhPR_Diam = 0.125;                     % [in], Exhaust Hose Diameter
ExhPR_Area = (ExhPR_Diam/2)^2*pi;       % [in^2], Exhaust Hose Cross-section Area
ExhPR_HDiam = ExhPR_Diam;               % [in], Hydraulic Diameter - ASSUMES CIRCULAR TUBING
% TODO: Pipe Friction and Heat Transfer, Override Values
%% Exhaust S Subsystem
% Reservoir Parameters
ExhRS_Pres = 1;                         % [atm], Reservoir Pressure
ExhRS_Temp = 273.15;                    % [K], Reservoir Temperature
ExhRS_Diam = 0.125;                     % [in], Exhaust Hose Inlet Diameter
ExhRS_Area = (ExhRS_Diam/2)^2*pi;       % [in^2], Exhaust Cross-section Area - UNUSED
% Pipe Parameters
ExhPS_Length = 1;                       % [in], Exhaust Pipe Length
ExhPS_Diam = 0.125;                     % [in], Exhaust Hose Diameter
ExhPS_Area = (ExhPS_Diam/2)^2*pi;       % [in^2], Exhaust Hose Cross-section Area - UNUSED
ExhPS_HDiam = ExhPS_Diam;               % [in], Hydraulic Diameter - ASSUMES CIRCULAR TUBING
% TODO: Pipe Friction and Heat Transfer, Override Values
%% Supply P Subsystem
% Reservoir Parameters (Need info on where compression is developed)
SupRP_Pres = 1;                         % [atm], Reservoir Pressure
SupRP_Temp = 273.15;                    % [K], Reservoir Temperature
SupRP_Diam = 0.25;                      % [in], Exhaust Hose Inlet Diameter
SupRP_Area = (SupRP_Diam/2)^2*pi;       % [in^2], Exhaust Cross-section Area
% Pressure Source Parameters
SupPS_Diff = 50;                        % [psi], Pressure differential
SupPS_DiamA = 0.25;                     % [in], Pressure source inlet diameter
SupPS_AreaA = (SupPS_DiamA/2)^2*pi;     % [in^2], Cross-section Area at inlet to pressure source
SupPS_DiamB = 0.25;                     % [in], Pressure source outlet diameter
SupPS_AreaB = (SupPS_DiamB/2)^2*pi;     % [in^2], Cross-section Area at outlet to pressure source
% Pipe Parameters
SupPP_Length = 1;                       % [in], Exhaust Pipe Length
SupPP_Diam = 0.25;                      % [in], Exhaust Hose Diameter
SupPP_Area = (SupPP_Diam/2)^2*pi;       % [in^2], Exhaust Hose Cross-section Area
SupPP_HDiam = SupPP_Diam;               % [in], Hydraulic Diameter - ASSUMES CIRCULAR TUBING
% TODO: Pipe Friction and Heat Transfer, Override Values
%% 5/2-Way Valve Subsystem
% Orifice AR Parameters
AR_MinDiam = 1e-5;                      % [in], Minimum restriction Diameter (Ideal = 0, but divide by zero)
AR_MinArea = (AR_MinDiam/2)^2*pi;       % [in^2], Minimum restriction Cross-section Area (Ideal = 0, but divide by zero)
AR_MaxArea = 16;                        % [mm^2], Maximum restriction Cross-section Area
AR_OriDiam = 0.25;                      % [in], Inlet and outlet Diameter
AR_OriArea = (AR_OriDiam/2)^2*pi;       % [in^2], Inlet and outlet Cross-section Area 
AR_Cd = 1;                              % Discharge coefficient (Ideal = 1)
% Orifice PA Parameters
PA_MinDiam = 1e-5;                      % [in], Minimum restriction Diameter (Ideal = 0)
PA_MinArea = (PA_MinDiam/2)^2*pi;       % [in^2], Minimum restriction Cross-section Area (Ideal = 0)
PA_MaxArea = 16;                        % [mm^2], Maximum resitrction Cross-section Area
PA_OriDiam = 0.25;                      % [in], Inlet and outlet Orifice Diameter 
PA_OriArea = (PA_OriDiam/2)^2*pi;       % [in^2], Inlet and outlet  Orifice Cross-section Area
PA_Cd = 1;                              % Discharge coefficient (Ideal = 1)
% Orifice PB Parameters
PB_MinDiam = 1e-5;                      % [in], Minimum restriction Diameter (Ideal = 0)
PB_MinArea = (PB_MinDiam/2)^2*pi;       % [in^2], Minimum restriction Cross-section Area (Ideal = 0)
PB_MaxArea = 16;                        % [mm^2], Maximum resitrction Cross-section Area
PB_OriDiam = 0.25;                      % [in], Inlet and outlet Orifice Diameter
PB_OriArea = (PB_OriDiam/2)^2*pi;       % [in^2], Inlet and outlet Orifice Cross-section Area
PB_Cd = 1;                              % Discharge coefficient (Ideal = 1)
% Orifice BS Parameters
BS_MinDiam = 1e-5;                      % [in], Minimum restriction Diameter (Ideal = 0)
BS_MinArea = (BS_MinDiam/2)^2*pi;       % [in^2], Minimum restriction Cross-section Area (Ideal = 0)
BS_MaxArea = 16;                        % [mm^2], Maximum resitrction Cross-section Area
BS_OriDiam = 0.25;                      % [in], Orifice Diameter, inlet and outlet
BS_OriArea = (BS_OriDiam/2)^2*pi;       % [in^2], Orifice Cross-section Area, inlet and outlet
BS_Cd = 1;                              % Discharge coefficient (Ideal = 1)
% Stroke-Area Conversion Parameters
% Assumes some interesting things about the spool:
%   -   That every hole is just barely not covered or just barely covered
%   -   Leak is equal for every hole
%   -   Conversion is equal for every hole
% TODO: Expand the model for individual orifice properties if these
% assumptions don't hold
% Needs stroke length, x leak offset, area leak offset, maximum stroke
% length
% Without this information, X_stroke MUST equal X_max (A binary ON/OFF
% model of the directional valve)
% X_Stroke is retrieved from the signal editor block
% X_Max therefore must be the maximum value in the signal editor and the
% absolute value of the maximum negative value
X_Max = 1;                              % [mm], Maximum spool stroke length
% X_Leak = ?;                           % [mm], Stroke leak offset
% A_Leak = ?;                           % [mm^2], Area leak offset
A_Max = 16;                             % [mm^2], Area maximum opening
SA_Gain = A_Max/X_Max;                  % [mm^2/mm], Distance converted to Area
% If all the required values are known, build the following formula in the
% Simulink
% [(X_Stroke - X_Leak)*(A_Max - A_Leak)/(X_Max - X_Leak)]+A_Leak
%% Presssure Pipes in the Main Model
% Pressure Pipe A Parameters
PA_Length = 12;                         % [in], Pipe Length
PA_Diam = 0.25;                         % [in], Pipe Diameter
PA_Area = (PA_Diam/2)^2*pi;             % [in^2], Pipe Cross-section Area
PA_HDiam = PA_Diam;                     % [in], Hydraulic Diameter - ASSUMES CIRCULAR PIPE
% Pressure Pipe B Parameters
PB_Length = 12;                         % [in], Pipe Length
PB_Diam = 0.25;                         % [in], Pipe Diameter
PB_Area = (PB_Diam/2)^2*pi;             % [in^2], Pipe Cross-section Area
PB_HDiam = PB_Diam;                     % [in], Hydraulic Diameter - ASSUMES CIRCULAR PIPE
%% Double-Action Actuator Model Subsystem
% Actuator Properties
ActuatorRod_Diam = 0.31;                            % [in], Rod Diameter
ActuatorRod_Area = (ActuatorRod_Diam/2)^2*pi;       % [in^2], Rod Cross-section Area
ActuatorPiston_Diam = 0.75;                         % [in], Piston Face Diameter
ActuatorPiston_Area = (ActuatorPiston_Diam/2)^2*pi; % [in^2], Rod Cross-section Area
ExtArea = ActuatorPiston_Area;                      % [in^2], Extension Pressure Interface
RetArea = ActuatorPiston_Area-ActuatorRod_Area;     % [in^2], Retraction Pressur Interface
ActuatorStroke_Length = 3.25;                       % [in], Actuator Rod Stroke length
% Chamber A Parameters
% Mechanical orientation: Negative
ChA_InterfaceArea = RetArea;                        % [in^2] Chamber A is the retraction chamber
ChA_InitialDisplacement = -ActuatorStroke_Length;   % [in] Piston starts fully retracted
ChA_DeadVol = 0.1*ChB_InterfaceArea*ActuatorStroke_Length;                                % [in^3] Dead volume in Chamber A (Ideal = 0, but divide by zero)
ChA_PortDiam = 0.1517;                              % [in] Port Diameter (From the minor diameter of 10-32 UNF)
ChA_PortArea = (ChA_PortDiam/2)^2*pi;               % [in^2] Port Area
% TODO: Get actual inlet minor diameter from hose insert
% Chamber B Parameters
% Mechanical orientation: Positive
ChB_InterfaceArea = ExtArea;                        % [in^2] Chamber B is the extension chamber
ChB_InitialDisplacement = ActuatorStroke_Length;    % [in] Piston starts fully retracted
ChB_DeadVol = 0.1*ChB_InterfaceArea*ActuatorStroke_Length;                                    % [in^3] Dead volume in Chamber B (Ideal = 0, but divide by zero)
ChB_PortDiam = 0.1517;                              % [in] Port Diameter (From the minor diameter of 10-32 UNF)
ChB_PortArea = (ChB_PortDiam/2)^2*pi;               % [in^2] Port Area
% TODO: Get actual inlet minor diameter from hose insert
% Mechanical Translational Blocks
% Damper, damping as a result of the piston's motion in the chamber
Actuator_Damp = 2000;                               % [N-s/m]
% Hard Stop
Actuator_StopU = ActuatorStroke_Length;             % [in] Hard stop Upper Bound
Actuator_StopL = 0;                                 % [in] Hard stop Lower Bound
% TODO: Develop reasonable values for contact stiffness and damping at the
% bounds
Actuator_RodMass = 1;                               % [kg] Actuator Rod Mass
% Load Subsystem
Const_Load = -0;                                   % [lbf] Load throughout the motion, negative opposes extension
% TODO: Make this load vary with position of the actuator rod
