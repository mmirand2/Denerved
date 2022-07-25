clear

rng(2); 

%function DriverBasic()

% Assign the patient data
data = Patient;

% Compute nominal parameter values and inital conditions
[gPars,Init] = load_global(data);
% pars = exp(gPars.pars); %Exponentiate the log-scaled nominal parameter values
pars = gPars.pars; %Exponentiate the log-scaled nominal parameter values

% Set up simulation time 
dt = 0.01;                                % Solution is displayed at time-step dt = 0.001
Hs= 60/pars(end);
T= round(1/Hs/dt)*dt;

gPars.period = 1;%period;                     % Cardiac cycle lenght added to data vector
data.Init     = Init;     
data.NumBeats =  300; %30                       % Number of cycles computed
% data.BeatsSS  =  70;                       % Results displayed after steady state has been obtained
tend          = data.NumBeats*gPars.period;         % Last time point
data.time     = [0:dt:tend];                % Time vector
data.Periods  = [0:gPars.period:data.NumBeats*gPars.period]; % Vector with times for periods
timeS         = [];
SolS          = [];
ptaS          = [];
eps  = 1e-6;

% t_R_STOP = 200;
% t_H_STOP = t_R_STOP+100;

% tup_initial = t_R_STOP;
% tend_initial = t_H_STOP;
% 
% td = t_H_STOP; %Don't tilt down
% ts   = 200;
% tu = tup_initial;
height = 25;
pars = [pars' height];

%%
k1    = 1;
k2    = round(T/dt)+k1;

while k2 < data.time(end)/dt %tend/dt+1+eps
    
    tdc = [k1-1:k2-1]*dt;
%   pars(end-2) = tup;

    options = odeset('AbsTol',gPars.ABS_TOL,'RelTol',gPars.REL_TOL);
    sol  = ode15s(@model_basic,[tdc(1) tdc(end)],data.Init,options,pars,tdc(1),T);
    %sol  = ode15s(@model_basic,[data.time(1) 10],data.Init,options,pars,gPars);
    sols = deval(sol,tdc);


    % Assigns each row as a temporary vector to store the solutions
    % For current time period
     timeS   = [timeS tdc(1:end-1)];
     SolS    = [SolS  sols(:,1:end-1)];

     data.Init    = SolS(:,end);

     %T       = round(1./sols(11,end)*(1+unifrnd(-1,1)*.02)/dt)*dt; %Adds noise
     T=1;
     k1      = k2; % Sets last index of this loop as first index for next loop
     k2 = k2+round(T/dt); % Redefining the last index
end
 
 %%
 plot( timeS,SolS')
 legend('V_rv','V_pa','V_pu','V_lv','V_tau','V_la','V_vl','V_tvu','Vtvl','Vtal','interpreter','none')
 hold on
 plot(t,x*10)
 grid

 plot(t,j)
 ylabel('Volume')
 xlabel('Time')
 legend('V_rv','V_pa','V_pu','V_lv','V_ta','V_la','V_vl','V_tv','H','interpreter','none')
%  plot(SolS([9,13],:)')

%%
% Generate time vecotr for displaying solution starting at t=data.BeatsSS (70s)
ID  = max(find(data.time    <= data.BeatsSS*gPars.period+1e-6));
IDp = max(find(data.Periods <= data.BeatsSS*gPars.period+1e-6));
timesol = data.time(ID:end);
timeper = data.Periods(IDp:end);

% Elastance function driver parameters
    H           = pars(1);                       % Elastance fctn param (1/s^2)
    T2          = pars(2);                       % Elastance fctn param (s)
 % Left ventricle free wall parameters
    E_lvf       = pars(3);                       % LV free wall elastance (kPa/mL)
    V_d_lvf     = pars(4);                       % LV end systolic zero pressure volume (mL)
    P_0_lvf     = pars(5);                       % LV end diastolic pressure (kPa)
    lambda_lvf  = pars(6);                       % LV end diastolic pressure (1/mL)
    V_0_lvf     = pars(7);                       % LV end diastolic pressure (mL)
    
    % Right ventricle free wall parameters
    E_rvf       = pars(8);                       % RV free wall elastance (kPa/mL)
    V_d_rvf     = pars(9);                       % RV end systolic zero pressure volume (mL)
    P_0_rvf     = pars(10);                      % RV end diastolic pressure (kPa)
    lambda_rvf  = pars(11);                      % RV end diastolic pressure (1/mL)
    V_0_rvf     = pars(12);                      % RV end diastolic pressure (mL)

      % Pulmonary artery and vein parameters
    E_pa        = pars(13);                      % PA elastance (kPa/mL)
    E_pu        = pars(14);                      % PU elastance (kPa/mL)
    R_pul       = pars(15);                      % Pulmonary vascular resistance (kPa*s/mL)
    P_th        = 0; %-pars(16);

     %Thoracic parameters
    E_ta        = pars(16);                      % TA elastance (kPa/ml)
    E_tv        = pars(17);                      % TV elastance (kPa/ml)
%   R_t         = pars(18);                      % Thoracic resistance (kPa*s/mL)

    % Lower body parameters
    E_la        = pars(18);                      % LA elastance (kPa/ml)
    E_vl        = pars(19);                      % LV elastance (kPa/ml)
%   R_l         = pars(21);                      % lOWER FLOW resistance (kPa*s/ml)
    R_lu        = pars(20);   
    R_ul        = pars(21);                      % UPPER-LOWER resistance (kPa*s/ml) 
    
    % Heart valve paramenters
    R_mt        = pars(22);                  % Mitral valve resist (kPa*s/mL)
    R_av        = pars(23);                      % Aortic valve resist (kPa*s/mL)
    R_tc        = pars(24);                      % Tricuspid vlv resist (kPa*s/mL)
    R_pv        = pars(25);                      % Pulmon vlv resist (kPa*s/mL)
        
    EdM         = pars(26);
    Edm         = pars(27);
    tauPm       = pars(28);
    tauE        = pars(29);
    kE          = pars(30);
    p2E         = pars(31);
    Es          = pars(32);

    HM          = pars(33);
    Hm          = pars(34);
    p2H         = pars(35);
    kH          = pars(36);
    tauH        = pars(37);
    Esr         = pars(38);
    EdMr        = pars(39);
    Edmr        = pars(40);
    p2Er        = pars(41);

    RlM         = pars(42);
    Rlm         = pars(43);
    p2Ra        = pars(44);
    RtM         = pars(45);
    Rtm         = pars(46);
    p2Ru        = pars(47);
    taur        = pars(48);
    kR          = pars(49);


%%% Ventricles 
V_lv        = SolS(4,:);                     % Left ventricular volume (mL)
V_rv        = SolS(1,:);                     % Right ventricular volume (mL)


%%% Left ventricular pressure
Ed  = SolS(9,:);
Es  = pars(32);

Ts = 0.001.*(0.82/1.82)*(522-1.87*60/T);
Tr = 0.001.*(1/1.82)*(522-1.87*60/T);
% 
e_t=[];
for i=1:length(Ed)
e_t(i) = ElastanceBasic(timeS(i)-ts,T,Ts,Tr,Ed(i),Es); % Ventriular Elastance function from Elastance1.m
end
P_lv= e_t.*V_lv;

%%% Right ventricular pressure
Edr  = SolS(13,:);
e_tr=[];
for i=1:length(Ed)
e_tr(i) = ElastanceBasic(timeS(i)-ts,T,Ts,Tr,Edr(i),Esr); % Ventriular Elastance function from Elastance1.m
end
P_rv= e_tr.*V_rv;
 

%%% Pulmonary arteries (PA)
V_pa        = sols(2,:);                    % Volume (mL)
P_pa        = E_pa * V_pa + P_th;           % Pressure (mmHg)

%%% Pulmonary veins (PU)
V_pu        = sols(3,:);                    % Volume (mL)
P_pu        = E_pu * V_pu + P_th;           % Pressure (mmHg)

%%% Systemic arteries (SA)
V_ta        = SolS(5,:);                    % Volume (mL)
P_ta        = E_ta * V_ta;                  % Pressure (mmHg)

%%%H
HR= SolS(11,:)*60;


%%% Thoracic veins (SV)
V_tv        = sols(8,:);                    % Volume (mL)
P_tv        = E_tv * V_tv;                  % Pressure (mmHg)

%%% Lower body veins
V_vl        = sols(7,:);                    % Volume (mL)
P_vl        = E_vl * V_vl;                  % Pressure (mmHg)

%%% Lowe body arterys
V_la        = sols(7,:);                    % Volume (mL)
P_la        = E_la * V_la;                  % Pressure (mmHg)



%%% Cardiac output (CO)
% Q_sys  = (P_lv - P_ta) / R_mt;             % Systemic cardiac output (mL/s)
Q_sys  = (P_ta - P_tv) / R_mt;   
Q_pul  = (P_pa - P_pu) / R_pul;             % Pulmonary cardiac output (mL/s)

%%% Information extracted over each cardiac cycle displayed (from t=70 to 80)
C1 = data.BeatsSS;
C2 = data.NumBeats;
PPC = gPars.pts;

Periods = data.Periods;
for i = C1:C2  
  I1 = (i-1)*PPC+1;
  I2 = i*PPC+1;
  tdc   = data.time(I1:I2);                 % Time points in current cycle
  P_sac = P_ta(I1:I2);                      % SA pressure in current cycle
  P_pac = P_pa(I1:I2);                      % PA pressure in current cycle
  P_rvc = P_rv(I1:I2);                      % RV pressure in current cycle
  V_lvc = V_lv(I1:I2);                      % LV volume in current cycle
  V_rvc = V_rv(I1:I2);                      % RV volume in current cycle
  
  PsaSys(i-C1+1) = max(P_sac);              % SA systolic pressure
  PsaDia(i-C1+1) = min(P_sac);              % SA diastolic pressure
  PpaSys(i-C1+1) = max(P_pac);              % PA systolic pressure
  PpaDia(i-C1+1) = min(P_pac);              % PA diastolic pressure
  PrvSys(i-C1+1) = max(P_rvc);              % Right ventriclular systolic pressure
  PrvDia(i-C1+1) = min(P_rvc);              % Right ventriclular diastolic pressure
  CO_sys(i-C1+1) = trapz(tdc,Q_sys(I1:I2))/(tdc(end)-tdc(1)); % Systemic cardiac output
  CO_pul(i-C1+1) = trapz(tdc,Q_pul(I1:I2))/(tdc(end)-tdc(1)); % Pulmonary cardiac output
  
  CO_Vl = (max(V_lvc) - min(V_lvc))/(tdc(end)-tdc(1)); % Cardiac output left ventricle
  CO_Vr = (max(V_rvc) - min(V_rvc))/(tdc(end)-tdc(1)); % Cardiac output right ventricle
  % Note all calculations of cardiac output should be the same.
  
  display([tdc(1) tdc(end) CO_sys(i-C1+1) CO_pul(i-C1+1) CO_Vl CO_Vr]) % Display values for current cardiac cycle
  Ppum(i-C1+1) = trapz(tdc,P_pu(I1:I2))/(tdc(end)-tdc(1));  % ??????
end;

% Data
P_SAsys   = data.P_SAsys *ones(size(PsaSys));   % Systolic aortic press (mmHg)    
P_SAdia   = data.P_SAdia*ones(size(PsaDia));   % Diastolic aortic press (mmHg)   
% P_PAsys   = data.P_PAsys *ones(size(PpaSys));   % Systolic PA press (mmHg)    
% P_PAdia   = data.P_PAdiast*ones(size(PpaDia));   % Diastolic PA press (mmHg)  
P_PAsys   = 21 *ones(size(PpaSys));   % Systolic PA press (mmHg)    
P_PAdia   = 8*ones(size(PpaDia));   % Diastolic PA press (mmHg)   
% P_RVsys   = data.P_RVsyst *ones(size(PrvSys));   % Systolic Right Ventricle press (mmHg)   
% P_RVdia   = data.P_RVdiast*ones(size(PrvDia));   % Diastolic Right Ventricle  press (mmHg)  
P_RVsys   = P_SAsys.*1.05.*ones(size(PrvSys));   % Systolic Right Ventricle press (mmHg)   
P_RVdia   = P_SAdia.*1.05.*ones(size(PrvDia));   % Diastolic Right Ventricle  press (mmHg)  
COd       = data.CO*ones(size(CO_sys));   % Cardiac output (mL/s)
% PCWp      = data.P_PCWave*ones(size(Ppum));      % Pulmonary capilary wedge pressure (mmHg)


% Total stressed volume (should be constant)
Vtotal = V_lv + V_rv + V_pa + V_pu + V_ta + V_tv + V_vl + V_la;

% Save results to file
% save Results233.mat;

% Figure(1) 
% Left ventricle, pulmonary vein, & aorta comparison of computed 
% results (solid lines) and data (broken lines) (2,3,1) 
%
% Right ventricle, pulmonary artery, & vena cava comparison of 
% computed results (solid lines) and data (broken lines) (2,3,2) 
%
% comparison of computed cardiac output and data. (2,3,3)
%
% comparison of left and right ventricular pressure. (2,3,4)
%
% computed left and right ventricular volume (no data available).(2,3,5)
%
% Left and right ventricular pressure volume loop (PV loop) with calculated 
% stroke work. (2,3,6)

t1 = timesol(1);
t2 = timesol(1)+5;
figure(1);
subplot(2,3,1);hold on;
    h=plot(timesol,P_lv(ID:end),timesol,P_ta(ID:end),timesol,P_pu(ID:end)); hold on;
    set(h(1),'Color',[.7 .06 .06])  
    set(h(2),'Color',[1 .50 .50])
    set(h(3),'Color',[1 .54 0])
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Pressure (mmHg)');
%     g=plot(timeper,PCWp,'--','linewidth',3);
%     set(g,'Color',[1 .54 0])
    g=plot(timeper,P_SAsys,'--',timeper,P_SAdia,'--','linewidth',3);
    set(g,'Color',[1 .50 .50])
    legend('P_{lv}','P_{ta}','P_{pv}');
    grid on;
    xlim([t1 t2]);    
    ylim([0 175]);
    
subplot(2,3,2);hold on;
    h = plot(timesol,P_rv(ID:end),'b',timesol,P_pa(ID:end),'c',timesol,P_tv(ID:end),'g'); hold on
    g = plot(timeper,P_PAsys,'c--',timeper,P_PAdia,'c--');hold on;
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Pressure (mmHg)');
    plot(timeper,P_RVsys,'b--',timeper,P_RVdia,'b--','linewidth',3);
    legend('P_{rv}','P_{pa}','P_{sv}');
    xlim([t1 t2]);  
    ylim([0 50]);
    grid on;
    
subplot(2,3,3);hold on;
    g=plot(timeper,COd,'--k',timeper,CO_sys,'--r',timeper,CO_pul,'--b','linewidth',3);
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Flow (mL/sec)');
    legend('CO_{data}','CO_{sys}','CO_{pul}');
    grid on;
    xlim([t1 t2]);
   
    
subplot(2,3,4);hold on;
    h=plot(timesol,P_lv(ID:end),timesol,P_rv(ID:end),'b');
    set(h(1),'Color',[.7 .06 .06])
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Pressure (mmHg)');
    legend('P_{lv}','P_{rv}');
    grid on;
    xlim([t1 t2]);
    ylim([0 165]);
    
subplot(2,3,5);hold on;
    h=plot(timesol,V_lv(ID:end),timesol,V_rv(ID:end),'b');
    set(h(1),'Color',[.7 .06 .06])
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Time (sec)');
    ylabel('Volume (mL)');
    legend('V_{lv}','V_{rv}');
    grid on;
    xlim([t1 t2]);
    ylim([40 165]);

subplot(2,3,6);hold on;
    h=plot(V_rv(ID:end),P_rv(ID:end),'b');
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    grid on;
    hold on;
    h=plot(V_lv(ID:end),P_lv(ID:end));
    set(h,'Color',[.7 .06 .06])
    set(h,'linewidth',3);
    set(gca,'fontsize',16);
    xlabel('Volume (mL)');
    ylabel('Pressure (mmHg)');
    grid on;
    legend('RV','LV')
    ylim([0 167]);
    
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run model 

global ODE_TOL
t_R_STOP = 200;
t_H_STOP = t_R_STOP+100;
dt   = .01;
tup_initial = t_R_STOP;
tend_initial = t_H_STOP;
td = t_H_STOP; %Don't tilt down
eps  = 1e-6;
ts   = 200;
height = 25;

kR = 23;
kH = 35;
BV = 3887.9;
[x0, Init] = load_global_clean(kH,kR,BV); 
tu = tup_initial;
x0 = [x0 tu td height];
pars = x0;



% Initialize empty vectors 
timeS = [];
TS    = []; 
pauS  = []; 
VvlS  = [];
pmS   = [];
HcS   = [];
RaupS = [];
RalpS = [];
EdS   = [];
CauS  = [];
CvuS  = [];
CvlS  = [];
SolS  = [];

H     = Init(end);
T     = round(1/H/dt)*dt;

k1    = 1;
k2    = round(T/dt)+k1;
tend = tend_initial;
tup  = tup_initial;
UP = 0;
while k2 < tend/dt+1+eps

    tdc = [k1-1:k2-1]*dt;
    pars(end-2) = tup;

    options = odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);
    sol     = ode15s(@model_clean,[tdc(1) tdc(end)],Init,options,pars,tdc(1),T); % Solve the ODEs
    sols    = deval(sol,tdc);

    % Assigns each row as a temporary vector to store the solutions
    % For current time period
     Vau     = sols(1,:);
     Vvl     = sols(4,:);
     pm      = sols(6,:);    
     Raup    = sols(7,:);
     Ralp    = sols(8,:);
     Ed      = sols(9,:);
     Cau     = pars(3);
     Hc      = sols(10,:);

     pau     = Vau./Cau;

     timeS   = [timeS tdc(1:end-1)];
     pauS    = [pauS  pau(1:end-1)];
     HcS     = [HcS   Hc(1:end-1)];
     SolS    = [SolS  sols(:,1:end-1)];

    Init    = sols(:,end);
%     T       = round(1./sols(10,end)/dt)*dt;
    T       = round(1./sols(10,end)*(1+unifrnd(-1,1)*.02)/dt)*dt; %Adds noise
    k1      = k2; % Sets last index of this loop as first index for next loop
    k2 = k2+round(T/dt); % Redefining the last index
end

  timeS     = [timeS tdc(end)]; 
  pauS      = [pauS  pau(end)];
  HcS       = [HcS   Hc(end)];
  SolS      = [SolS  sols(:,end)];
 

%% Plotting
fontsize = 15;
co = 'k';

ylh = [.85,1.15];
ylp = [50,130];

tbefore = 75;
tafter = 75;
tstartS = tup-tbefore;
tendS = tup + tafter;
t_tilt = tup;%750;
tilt_indm = find(abs(timeS-t_tilt) == min(abs(timeS-t_tilt)));%Tilt ind
ind1 = find(abs(timeS-tstartS) == min(abs(timeS-tstartS))); %Start ind
ind2 = length(timeS);%find(abs(timeS-tendS) == min(abs(timeS-tendS))); %End ind
sm1 = ind1:ind2;%ind span for simulation


figure(1)
clf
subplot(2,1,1)
hold on
plot(timeS(sm1),HcS(sm1),'r','linewidth',2) %Plot model sim
plot(ones(2,1).*timeS(tilt_indm),ylh,'--k','linewidth',2)  %Tilt line
xticks([])
ylabel('H (bps)')
set(gca,'fontsize',fontsize)
xlim([tstartS,tendS])
% ylim(ylh)

subplot(2,1,2)
hold on
plot(ones(2,1).*timeS(tilt_indm),ylp,'--k','linewidth',2)
plot(timeS(sm1),pauS(sm1),'r')
plot(timeS(sm1),SolS(6,sm1),'k','linewidth',2)
ylabel('P_{au} (mmHg)')
xlabel('Time (s)')
xlim([tstartS,tendS])
set(gca,'fontsize',fontsize)








