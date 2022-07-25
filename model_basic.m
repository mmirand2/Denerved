function xdot= model_basic(time,y,pars,ts,T)

   % Unpack the state variable vector
    V_rv  = y(1); %Rigth ventricule  
    V_pa  = y(2); %Pulmonary arteries
    V_pv  = y(3); %Pulmonary veins
    V_lv  = y(4); %Left ventricule
    V_tau = y(5); %Thoracic arteries
    V_la  = y(6); %lower body arteries
    V_vl  = y(7); %Lower body veins
    V_tvu = y(8); %Thoracic Veins
    V_tvl = y(9);
    V_tal = y(10);
    


    %% Parameters

      % Pulmonary artery and vein parameters
    E_pa        = pars(1);                      % PA elastance (kPa/mL)
    E_pu        = pars(2);                      % PU elastance (kPa/mL)
    R_pul       = pars(3);                      % Pulmonary vascular resistance (kPa*s/mL)
  
     %Thoracic parameters
    E_tau        = pars(4);                      % TA elastance (kPa/ml)
    E_tvu        = pars(5);                      % TV elastance (kPa/ml)

    % Lower body parameters
    E_la        = pars(6);                      % LA elastance (kPa/ml)
    E_vl        = pars(7);                      % LV elastance (kPa/ml)
    R_l         = pars(8);                      % lOWER FLOW resistance (kPa*s/ml)
    R_lu        = pars(9);   
    R_ul        = pars(10);                      % UPPER-LOWER resistance (kPa*s/ml) 
    R_tu        = pars(11);

    % Heart valve paramenters
    R_mt        = pars(12);                  % Mitral valve resist (kPa*s/mL)
    R_av        = pars(13);                      % Aortic valve resist (kPa*s/mL)
    R_tc        = pars(14);                      % Tricuspid vlv resist (kPa*s/mL)
    R_pv        = pars(15);                      % Pulmon vlv resist (kPa*s/mL)
        
    Es          = pars(16);
    Esr         = pars(17);

    R_tul       = pars(18);
    E_tvl       = pars(19);
    E_tal       = pars(20);
    R_tlu       = pars(21);
    VM_vl       = pars(22);
    m_vl        = pars(23);

    R_tl        = pars(24);
    A_u         = pars(25);
    B_u         = pars(26);
    A_l         = pars(27);
    B_l         = pars(28);
    Tinsp       = pars(29);
    Tc          = pars(30);
    H           = pars(31);
    Ed          = pars(32);
    Edr         = pars(33);

   % respiration
    
    tresp= mod(time,Tc);

    if tresp< Tinsp
        P_extu= (A_u/2)*( cos((2*pi*tresp)/Tinsp)-1 )+B_u;  %thorax
    else
        P_extu= B_u;
    end

    if tresp< Tinsp
        P_extl= (A_l/2)*( cos((2*pi*tresp)/Tinsp)-1 )+B_l;  %abdomen
    else
        P_extl= B_l;
    end

    
    % Elastance function left
    % Heart parameters

    P_tau = E_tau * (V_tau);                                                %Blood Pressure: Upper Thorax 
    Ts    = 0.001.*(0.82/1.82)*(522-1.87*60/T);                                %Time start sistole
    Tr    = 0.001.*(1/1.82)*(522-1.87*60/T);                                   %Time start diastole
    
    %carotid pressure
   
    e_t = ElastanceBasic(time-ts,T,Ts,Tr,Ed,Es); % Ventriular Elastance function from Elastance1.m
    P_lv= (e_t)*V_lv;
   
   % Elastance function rigth
   % Heart parameters

    e_tr = ElastanceBasic(time-ts,T,Ts,Tr,Edr,Esr); % Ventriular Elastance function from Elastance1.m
    P_rv= (e_tr)*V_rv; %0.55
    
    % Pulmonary pressure and flow
    P_pa = E_pa * (V_pa);     % PA pressure (mmHg)
    P_pu = E_pu * (V_pv);     % PU pressure (mmHg)
    Q_pul = (P_pa-P_pu) / R_pul;     % Pulmonary flow (mL/s)

    %%
    % Thoracic pressure and flow upper
    P_tal   = E_tal * (V_tal); 
    P_tvl   = E_tvl * (V_tvl); 

    
    P_tvu   = E_tvu * (V_tvu);             % Toracic veins pressures (mmHg) ***Agregar modulacion respiratoria
%     P_ta  = E_ta * (V_ta);
    Q_tu    = (P_tau-P_tvu)/R_tu;          % Thoracic Flow upper           (ml/s) ***Aca mirar info sobre el intercambio en torax
    Q_tul   = (P_tau-P_tal)/R_tul;
    Q_tlu   = (P_tvl-P_tvu)/R_tlu;
    

    % Thoracic pressure and flow lower
     Q_tl   =(P_tal-P_tvl)/R_tl; %reemplazar luego con R_tl  


    %%
    % Lower Body aretery pressure and flow

    P_la= E_la * V_la;
    P_vl= (1/m_vl) * log(VM_vl/(VM_vl - V_vl));
    
    Q_ul = (P_tal-P_la)/R_ul;
    Q_l  = (P_la-P_vl)/R_l;

    % Lower body veins flow
    

    if P_vl>P_tvl
    Q_lu = (P_vl-P_tvl)/R_lu;
    else
        Q_lu =0;
    end


%%

    
    % Valve flow
    if P_pu > P_lv                      % Flow through mitral valve (mL/s)
        Q_mt = (P_pu - P_lv)/(R_mt);
    else
        Q_mt = 0;
    end
    if P_lv > P_tau                      % Flow through aortic valve (mL/s)
        Q_av = (P_lv - P_tau)/R_av;
    else
        Q_av = 0;
    end
    if P_tvu > P_rv                      % Flow through tricuspid valve (mL/s)
        Q_tc = (P_tvu-P_rv)/R_tc;
    else
        Q_tc = 0;
    end
    if P_rv > P_pa                      % Flow through pulmonary valve (mL/s)
        Q_pv = (P_rv - P_pa)/R_pv;
    else
        Q_pv = 0;
    end
%% Changues in volume
dVrv  = Q_tc-Q_pv;            % Volume right ventricule (tricuspide-pulmonary valve)
dVpa  = Q_pv-Q_pul;           % Volume pulmonary artery (pulmonary valve- pulmonary flow)
dVpv  = Q_pul-Q_mt;           % Volume pumonary veins (pulmonary flow-mitral valve)
dVlv  = Q_mt-Q_av;            % Volume left heart (mitral-aortic valve)

dVtau = Q_av-Q_tul-Q_tu;     % Volume thorax artery upper*
dVtal = Q_tul-Q_tl-Q_ul;     % Volume thorax artery lower (LF- thorax flow-uppertolower) *

dVtvl = Q_lu+Q_tl-Q_tlu;     % Volume thorax vein lowe  (lower to upper+thoracic-systemic vein)*
dVtvu = Q_tlu+Q_tu-Q_tc;    % Volume thorax vein upper

dVla  = Q_ul-Q_l;             % Volume lower artery (uppertolower-lower flow)
dVvl  = Q_l-Q_lu;             % Volume lower vein (lower flow- lowertoupper)



xdot=[dVrv,dVpa,dVpv,dVlv,dVtau,dVla,dVvl,dVtvu,dVtvl,dVtal]';



%Regional Distribution of Cardiac Output: Normal Values in Man Determined
%by Video Dilution technique: torax diafragma

%Anatomy, Thorax, Internal Mammary (Internal Thoracic) Arteries

%Anatomy, Abdomen and Pelvis, Arteries and Veins
































