function [gPars, Init] = load_global(data)

BW        = data.BW;                         % Body weight (kg)                
Hgt       = data.Hgt;                        % Height (cm)                     
Gender    = data.Gender;                     % Gender (1=female, 2=male)       

%% Total blood volume (mL)
if (Gender == 2)  
        TotBV = ((0.3669 * (Hgt/100)^3) + (0.03219 * BW) + 0.6041) * 1000; % Total blood volume (mL)
        BSA   = 0.000579479 * BW^0.38 * Hgt^1.24;                          % Body surface area (m^2)
    else
        TotBV = ((0.3561 * (Hgt/100)^3) + (0.03308 * BW) + 0.1833) * 1000; % Total blood volume (mL)
        BSA   = 0.000975482 * BW^0.46 * Hgt^1.08;                          % Body surface area (m^2)
end

%% Cardiac Output (mL/s)

Ave_HR        = data.Ave_HR;                     % Average heart rate (beats/min)
gPars.period  = 60/Ave_HR;                   % Period of heart beat (s)
CO            = data.CO;                  % Cardiac output (thermodilution) (mL/s) Finapres
TotFlow       = CO;

HI            = 60/Ave_HR;
%% Stroke Volume
SV    = CO*(gPars.period);   % Left stroke volume (mL/beat) % (mL/s * s/beat) = mL/beat)

%% Tolerance for ODE solver
gPars.ABS_TOL  = 1e-8;
gPars.REL_TOL  = 1e-8;

%% Blood Pressure

% Systemic Arteries (SA)
P_TAsys  = data.P_SAsys;                   % Systolic SA pressure  (mmHg)    
P_TAdia  = data.P_SAdia;                  % Diastolic SA pressure (mmHg)   
P_TAm    = P_TAsys/3 + 2*P_TAdia/3;      % SA mean pressure      (mmHg)

P_TALsys = P_TAsys*0.98;
P_TALdia = P_TAdia*0.98;
P_TALm   = P_TALsys/3 + 2*P_TALdia/3;

P_LAsys  = data.P_SAsys*0.95;                   % Systolic SA pressure  (mmHg)    
P_LAdia  = data.P_SAdia*0.95;                  % Diastolic SA pressure (mmHg)   
P_LAm    = P_LAsys/3 + 2*P_LAdia/3;      % SA mean pressure      (mmHg)

% Systemic Veins (SV)
P_TVU      = 6;
P_TVL      = 8; %8*****                   % Diastolic SV pressure (mmHg) Understanding basic vein physiology and venous blood pressure through simple physical assessments
P_VL       = 10;

% Pulmonary Arteries (PA)
P_PAsys   = 24;%21;                         % Systolic  PA pressure (mmHg)  
P_PAdia   = 8;                              % Diastolic PA pressure (mmHg)  

P_PAmean   =P_PAsys/3 + 2*P_PAdia/3;

% Pulmonary Veins (PU)
P_PV       = 5;                              % Pulmonary vein pressure (mmHg) 

% Left ventricle (LV)
% P_LV       = 2.5;  %      
P_LVsyst  =  1.01 * P_TAsys;                 % Systolic  LV pressure (mmHg)
P_LVdia   =  1.01*P_TAdia;


% Right Ventricle (RV)
P_RVsyst   = P_PAsys*1.05;                   % Systolic  RV pressure (mmHg) BC Lampert,  
P_RVdia    = P_PAdia*1.05;                      %Recordar que esto es por ahora
                            

V_LVM      = 78*BSA;%93*BSA-16;%78*BSA;        % Max LV volume (End Diastolic Volume in LV (mL) 
V_RVM      =(78*BSA)*0.95; %0.9*V_LVM;%78*BSA;        % Max RV volume (10% higher than the  V_LV)  


V_LVm      = V_LVM - SV;       % Min LV volume (End Systolic Volume in LV mL)
V_RVm      = V_RVM-SV;        % Min RV volume (10% higher than the  V_LV)

V_d_lvf    = 10;                    % 4 LV end systolic zero pressure volume (mL) from Smith & Andreassen
V_d_rvf    = 0.9*V_d_lvf;           % 9 RV end systolic zero pressure volume (mL) from Smith & Andreassen


% Distribution of volume outside the heart (fractions add to 1)
sart = 0.13;     % Benekin 
part = 0.03;     % Benekin 
svein= 0.65;     % Benekin  (note we took volume from the systemic veins and put into the heart)
pvein= 0.11;     % Benekin  

% Circulating (stressed volume)
Circ_pa =  0.58 * part*TotBV;   % PA volume from Benekin   
Circ_pu =  0.11 * pvein*TotBV;  % PU volume from Benekin  
Circ_sa =  0.27 * sart*TotBV;   % SA volume from Benekin
Circ_sv =  0.08 * svein*TotBV;  % SV volume from Benekin

% Total circulating (stressed volume)
CircBV =  Circ_pa + Circ_pu + Circ_sa + Circ_sv;


%%
% Flows
qtu  = (TotFlow*0.80)*0.35; % Upper body flow, arteries --> veins (Upper)
qtul  = (TotFlow*0.80)*0.65;
qul  = TotFlow*0.20;        % Upper body arteries --> lower body arteries
qlu  = qul;                 % Lower body veins --> upper body veins
ql   = qul;                 % Lower body flow, arteries --> veins
qtl  = qtul;
qtlu = qtul;

%% Resistences

% Peripheral resistances (Ohm's law)
R_tu  = (P_TAdia-P_TVU)/qtu;  % Thoracic resistance
R_tul = (P_TAm-P_TALm)/qtul;
R_tl  = (P_TALm-P_TVL)/qtl;%(P_TALdia-P_TVU)/qtl;
R_l   = (P_LAdia-P_VL)/ql;  % lower body resistance 
R_ul  = (P_TAm-P_LAm)/qul;
R_lu  = (P_VL-P_TVL)/qlu;
R_tlu = (P_TVL-P_TVU)/qtlu;

% Valve resistances
R_mt  = 0.0025;%(P_PVdia-P_LVdia)/CO;%0.0025;%0.25./CO;        % 19 Mitral valve resistance (mmHg*s/mL)
R_av  = (P_LVsyst  - P_TAsys)/CO;         % 20 Aortic valve resistance (mmHg*s/mL)
R_tc  = 0.0025; %(P_TVdia-P_RVdia)/CO; %0.0011;%0.25./CO;        % 21 Tricuspid valve resistance (mmHg*s/mL)
R_pv  = (P_RVsyst  - P_PAsys)/CO;         % 22 Pulmonary valve resistance (mmHg*s/mL)
     
% Pumonary resistances
R_pul   = (P_PAdia - P_PV)/CO;           % 15 Pulmonary vascular resistance (mmHg*s/mL) 

%% Elastances

% Pulmonary artery and vein parameters
E_pa = P_PAsys/Circ_pa;                % 13 PA artery elastance (mmHg/mL)  
E_pu = P_PAmean/Circ_pu;                   % 14 PU elastance (mmHg/mL)  

% Systemic artery and vein parameters
E_tau    = P_TAsys/(Circ_sa*0.80)*0.35;                % 16 SA elastance (mmHg/mL)
E_tvu    = P_TVU/(Circ_sv*0.80)*0.35; 

E_la    = P_LAsys/(Circ_sa*0.20);
E_vl    = P_VL/(Circ_sv*0.20); %REVISAR SI ES IGUAL A Pvl
                                                % 17 SV elastance (mmHg/mL)
E_tal   = P_TALsys/(Circ_sa*0.80)*0.65;
E_tvl   = P_TVL/(Circ_sv*0.80)*0.65; 



% Left ventricle free wall parameters


%%NUEVO

Es     = P_TAsys/V_LVm;                     % Maximum elastance of the heart [ EdM = part / (min Vlv - Vd)
Ed     = P_TAdia/V_LVM;

%% elastancia rigth

Esr     = P_PAsys/V_RVm;
Edr     = P_TALdia/V_RVM;

%% Respiration
A_u= 0;
B_u=  0;

B_l= 0;
A_l=  -1.5;

Tinsp= 1.5;
Tc=8/3*(Tinsp);

%%
VMvl = 4*(Circ_sv*0.20);
m_vl  = log(VMvl/(VMvl - Circ_sv*0.20))/P_VL  ;  %Cvl/(VM_vl-Vvl)


%HR


   
%%%%% Initial conditions (stressed volumes) %%%%%
V_rv0 = V_RVM - V_d_rvf;                  % RV volume (mL)
V_pa0 = Circ_pa;                          % PA volume (mL)
V_pu0 = Circ_pu;                          % PU volume (mL)
V_lv0 = V_LVM - V_d_lvf;                  % LV volume (mL)
V_tau0 = Circ_sv*0.80*0.35;               % TV volume (ml)
V_la0  = Circ_sa*0.20;                    % LA volume (ml)  
V_vl0 = Circ_sv*0.20;                     % VL volume (ml)
V_tvu0 = Circ_sv*0.80*0.35;               % TV volume (ml)
V_tvl0 = Circ_sa*0.80*0.65;               % TA volume (ml)
V_tal0 = Circ_sa*0.80*0.65;               % TA volume (ml)

% V_sa0 = Circ_sa;                          % SA volume (mL)
% V_sv0 = Circ_sv;                          % SV volume (mL)

% Init = [V_lv0, V_rv0, V_pa0, V_pu0,V_tv0,V_vl0,V_ta0,V_la0,EdI,pm,HI,ppm,EdIr,R_l,R_tu];
Init = [V_rv0,V_pa0, V_pu0,V_lv0,V_tau0,V_la0,V_vl0,V_tvu0,V_tvl0,V_tal0];
    
%%%%% Model parameters %%%%%
pars = [                                 ...  
E_pa  E_pu R_pul                         ...  % 1 - 3
E_tau E_tvu                              ...   % 4 - 5
E_la E_vl R_l                            ...     % 6 - 7
R_lu R_ul R_tu R_mt R_av R_tc R_pv       ...     % 8 - 13 
Es                                       ...     % 14- 20
Esr                                      ...     % 21- 29
R_tul E_tvl E_tal R_tlu VMvl m_vl        ...     % 38- 43
R_tl                                     ...     % 44   
A_u B_u A_l B_l Tinsp Tc HI              ...
Ed Edr]';                   % 45-50             


%%%%% Structure with log-scaled model parameters
gPars.pars = pars;


end 