% This file contains the differential equations of the model and also
% calculates the new left ventricular pressure by calling the function
% Elastance.

% Big change: Cvu now a diff eq

function xdot = model_clean(t,y,pars,ts,T)

global Rav Rmv

% Define variables
Vau    = y(1);
Vvu    = y(2);
Val    = y(3); 
Vvl    = y(4); 
Vlv    = y(5);
pcm    = y(6);
Raup   = y(7);
Ralp   = y(8);
Ed     = y(9);
Hc     = y(10);

% Define resistances
Ral   = pars(1);
Rvl   = pars(2);    

% Define compliances
Cau = pars(3);
Cal	= pars(4);
Cvu = pars(5);

% Pressures
pau = Vau./Cau;
pal = Val./Cal;

VM_vl = pars(29);
m_vl = pars(30);
pvl = (1/m_vl) * log(VM_vl/(VM_vl - Vvl));
pvu = Vvu./Cvu;

% Heart parameters
Ts = 0.001.*(0.82/1.82)*(522-1.87*60/T);                                    
Tr = 0.001.*(1/1.82)*(522-1.87*60/T);                                       
Es     = pars(8);
Vd     = pars(9);

% Tilt Function
rhogh  = tiltftn_timescaled(t,pars);

%Carotid pressure
tilde_pars = pars;
tilde_pars(end) = 20; %Distance from aorta and carotid
rhogh  = tiltftn_timescaled(t,tilde_pars);
pc = pau - rhogh;

Elv  = ElastanceBasic(t-ts,T,Ts,Tr,Ed,Es); % Ventriular Elastance function from Elastance1.m
plv  = Elv*Vlv; % Stressed ventricular volume & Vd(consant), ventricular volume at zero dyastolic pressure

%Calculate flows
if plv > pau
   qav  = (plv - pau)/Rav;
else
   qav = 0;
end
if pvu > plv
  qmv  = (pvu - plv)/Rmv;
else
  qmv = 0;
end
qal  = (pau - pal + rhogh)/Ral;
qaup = (pau - pvu)/Raup;
qalp = (pal - pvl)/Ralp;
if pvl > pvu
    qvl  = (pvl - pvu - rhogh)/Rvl;
else
    qvl = 0;
end

% Mean pressure
tauPm = pars(31);
dpcm  = (pc-pcm)/tauPm;
% tilde_pars(end) = 20; %htilde is 20 - distance from aorta to carotid
% rhogh_tilde = tiltftn_timescaled(t,tilde_pars);
% dpm  = ((pau-rhogh_tilde)-pm)/tauPm;
% B = 0.5; %Respective contributions of carotid and aortic receptors
% dpm = (pau-B*rhogh_tilde-pm)/tauPm;

% Define resistances
kR   = pars(10);
taur = pars(11);

RaupM= pars(12);
Raupm= pars(13);
p2Ru = pars(14);

RalpM= pars(15);
Ralpm= pars(16);
p2Ra = pars(17);

% Ralpf = (RalpM - Ralpm)*p2Ra^kR/(pm^kR + p2Ra^kR) +Ralpm;
Ralpf = (RalpM - Ralpm)*p2Ra^kR/(pcm^kR + p2Ra^kR) +Ralpm;
dRalp = (-Ralp+Ralpf)/taur;

% Raupf = (RaupM - Raupm)*p2Ru^kR/(pm^kR + p2Ru^kR) +Raupm;
Raupf = (RaupM - Raupm)*p2Ru^kR/(pcm^kR + p2Ru^kR) +Raupm;
dRaup = (-Raup+Raupf)/taur;

%Cardiac contractility (decreasing minimum elastance)
kE  = pars(18);
tauE= pars(19);
EdM = pars(20);
Edm = pars(21);
p2E = pars(22);

% Edf = (EdM - Edm)*pm.^kE./(pm.^kE + p2E^kE) +Edm;
Edf = (EdM - Edm)*pcm.^kE./(pcm.^kE + p2E^kE) +Edm;
dEd = (-Ed + Edf)/tauE;

%heart rate
kH  = pars(23);
tauH= pars(24);
HM  = pars(25);
Hm  = pars(26);
p2H = pars(27);

% Hf = (HM-(Hm))*p2H^kH./(pm.^kH+p2H^kH) + Hm;
Hf = (HM-(Hm))*p2H^kH./(pcm.^kH+p2H^kH) + Hm;
dHc = (-Hc + Hf)/tauH;



xdot = [qav - qal - qaup; ...   %dVau upper arteries
        qvl + qaup - qmv; ...   %dVvu upper veins
        qal - qalp; ...         %dVal lower arteries
        qalp - qvl; ...         %dVvl lower veins
        qmv - qav ; ...         %dVlv left heart
        dpcm; dRaup; dRalp; dEd; dHc];%dCvl];
 
% xdot(end-3:end) = zeros(length(xdot(end-3:end)),1); %Test no control
% xdot(end-3:end-1) = zeros(length(xdot(end-5:end-1)),1); %Test no control
