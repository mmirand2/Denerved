% This function initializes the parameters for the model and sets initial
% values for the variables.R_t

%--------------------------------
%Changes from original:
% Ed -> Em
% Es -> EM
%--------------------------------

function [x0, Init] = load_global_clean(kH,kR,BV)

global ODE_TOL
global Rav Rmv

ODE_TOL  = 1e-8;

delta = 2.5;
%h = 167;
%BMI = 22;
%BV = 0.4948*BMI^0.425*h^1.575 - 1954;
TotalVol = BV;
TotFlow = TotalVol/60;     % Assume blood circulates in one min Cardiac output [ml/sec]
COd     = TotFlow*60/1000; % Steady state only [ l/min ]

HI      = .96;  % initial heart rate [sec]

% Flows (related to subject) 
qaup = TotFlow*0.80; % Upper body flow, arteries --> veins
qal  = TotFlow*0.20; % Upper body arteries --> lower body arteries
qvl  = qal;          % Lower body veins --> upper body veins
qalp = qal;          % Lower body flow, arteries --> veins

% Pressures (related to subject)     
pauD   = 80;  %From Data diastolic pressure
pauS   = 120; %From Data systolic pressure
pau    = (2/3)*pauD+(1/3)*pauS;  % Mean systolic pressure data (Upper body arteries)
pm     = pau;        % Initial mean pressure - Will oscilate around this value.

pal    = pau*0.99;   % Lower body arteries
pvu    = 2.75;       % Upper body veins - Note higher venous pressure. Should reflect pressures in pulmonary veins!
pvl    = 3.00;       % Lower body veins
plvD   = 2.5;        % Left ventricle diastolic
Vd     = 10;         % Left ventricular unstressed (residual) volume
Vlvm    = 50-Vd;      % Minimum left ventricular volume minus dead volume 50 average diastolic left ventricular volume
VlvM   = 110-Vd;     % Minimum left ventricular volume minus dead volume 50 average diastolic left ventricular volume

% Resistances (Ohm's law)
RaupI = (pau-pvu)/qaup;  % Upper body peripheral resistance
Ral   = (pau-pal)/qal;   % Upper body arteries --> lower body arteries
Rvl   = (pvl-pvu)/qvl;   % Lower body veins --> upper body veins
RalpI = (pal-pvl)/qalp;  % Lower body peripheral resistance
Rav   = 0.0001;          % Valve resistance
Rmv   = 0.0001;

TotalVolSV = TotalVol*0.85; % 85 percent in systemic circulation
Vart = TotalVolSV*.15; % 15 percent arteries and arterioles
Vven = TotalVolSV*.85; % 85 percent veins, capillaries, central veins

% Volume distribution (Beneken and deWit)
Vau = Vart*0.80*.30;     % 80 percent in upper body arteries of these 70 percent is unstressed
Val = Vart*0.20*.30;     % 20 percent in lower body arteries of these 70 percent is unstressed
Vvu = Vven*0.80*.075;    % 80 percent in upper body veins of these 92.5 percent is unstressed
Vvl = Vven*0.20*.075;    % 20 percent in lower body veins of these 92.5 percent is unstressed
                  
% Compliances, stressed volume percentages from Beneken are weighted
Cau = Vau/pauD;%pau;%120;%Vau/pau;  % Upper body artery compliance 
Cvu = Vvu/pvu;  % Upper body venous compliance
Cal  = Val/pal;  % Lower body artery compliance 
Cvl  = Vvl/pvl;  % Lower body venous compliance

% Ventricle parameters 
Tsfrac = 0.12;     % Time fraction for systolic phase, increasing elasticity 
Trfrac = 0.14;     % Time fraction for Systolic phase, decreasing elasticity
% % % EmI    = 0.0273;   % Minimum elastance of the heart [ Em = pveins / (Max Vlv - Vd) = 3 / (120-10) ]
% % % EM     = 3.6667;   % Maximum elastance of the heart [ EM = part / (min Vlv - Vd) = 120/(40-10) ]  % 120 is max arterial pressure
EdI     = plvD/VlvM;  % Minimum elastance of the heart [ Ed = pveins / (Max Vlv - Vd) 
Es     = pauS/Vlvm;   % Maximum elastance of the heart [ EdM = part / (min Vlv - Vd)

VauI = Vau;
VvuI = Vvu;
ValI = Val;
VvlI = Vvl;
VlvI = VlvM;%120-Vd; % Initial left ventricular volume is Vlv,max
               % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2657902/#!po=34.8485
               
%Resistance Upper Body Pars
taur   = delta*5; 
% kR     = 7; %Currently an input

RaupM = 3*RaupI; %Cite Miller for this, this is conservative in comparison I think radius -> 0.75 & -> 1.5, what rep?
Raupm = 0.2*RaupI; %With time scale, use 1.2 & .6 for time scaling 
alpha = (RaupI-Raupm)/(RaupM-Raupm);
p2Ru  = (pm^kR * alpha/(1 - alpha))^(1/kR); 

%Resistance Lower Body Pars
RalpM  = 3*RalpI;
Ralpm  = 0.2*RalpI; 
alpha  = (RalpI-Ralpm)/(RalpM-Ralpm);
p2Ra   = (pm^kR * alpha/(1 - alpha))^(1/kR);

%Contractility (Em)
EdM  = 1.25*EdI; 
Edm  = 0.01*EdI; 
tauE = delta*5;   
kE = 7; %20? Make higher for pulse pressure
alpha= (EdI-Edm)/(EdM-Edm);
p2E  = (pm^kE*(1 - alpha)/alpha)^(1/kE); 

% Lover body volumes
VMvl = 4*Vvl; %Choosen based on BP drop without controls 
mvl  = log(VMvl/(VMvl - Vvl))/pvl;  %Cvl/(VM_vl-Vvl)

%Cardiac cycle [Beregne HI fra p2H givet]
HM   = 200/60; %Hardcode to 200 bpm (Nes)
Hm   = .3;% Hardcode to 18 bpm (estimate)
% kH = 10;%Currently an input
tauH = delta*2.5;
alpha= (HI-Hm)/(HM-Hm);
p2H  =(pm^kH * alpha/(1 - alpha))^(1/kH);

Init = [VauI VvuI ValI VvlI VlvI pm RaupI RalpI EdI HI];

tauPm = delta*1; %NEW TO THIS SCALED TIME

% Parameter vector
x0 = [Ral Rvl ...                %1-2
      Cau Cal Cvu...       %3-5
      Tsfrac Trfrac Es Vd ...    %6-9
      kR taur RaupM Raupm p2Ru RalpM Ralpm p2Ra ... %10-17 
      kE tauE EdM Edm p2E...     %18-22  
      kH tauH HM Hm p2H,TotalVol,... %23-28 Added TotalVol
      VMvl, mvl,tauPm]; %29-31
  
  
