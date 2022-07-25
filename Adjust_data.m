%% Control
dat = load('C64549_fig_dat.txt');

Fs = round(1/mean(diff(dat(:,1)))); % 1/0.001 = 100
t = dat(:,1);
hr = dat(:,2);
p = dat(:,3);
tilt = 300; %Time in data when tilt happens
tilt_length = 7; %Taken from data
tilt_l_ind = 10*Fs;
ind_tilt = find(abs(t-tilt) == min(abs(t-tilt)),1);

rho       = 1.06;   % Density of blood in g/cm^3
g         = 982;    % Gravitational acceleration in cm/s^2
                    % Smithsonian physical tables, Smithsonian Institution,
                    % Seventh Revised Edition, prepared by Frerick E Fowle,
                    % Vol 71, #4, 1921


h  =  20; 
theta     = 60;             % Maximal tilt-angle, in degrees.
a         = 60/tilt_length;           % Tilt-speed
conv      = 1333.22;        % conversion cgs (cm,grams,sec) units to mmHg
arg = zeros(size(t));  
for z = 1:length(t)
    if t(z) < tilt
       arg(z) = 0;
    elseif t(z) < tilt+tilt_length
       arg(z) = a.*(t(z)-tilt);
    else
        arg(z) = theta; %No tilting down
    end
end

% Calculate the hydrostatic pressure in mmHg, 1Pa=0.0075mmHg
rhogh = rho.*g.*h.*(sin(arg.*pi./180)./conv); 

%

new_pressure = (p-80)*.8+90;%Line up the pulse pressures and subtract hydrostatic pressure
new_hr = (hr./60-.12); %Line up heartrates

writematrix([t,new_hr,new_pressure],'Adjusted_control_pp_hr.txt')



%% POTS
% clear
dat = load('20253_fig_dat.txt');

t = dat(:,1);
Fs = 1/mean(diff(dat(:,1))); 
p = dat(:,4);
tilt = 1840; %Time in data when tilt happens
tilt_length = 7; %Taken from data
tilt_l_ind = 10*Fs;
ind_tilt = find(abs(t-tilt) == min(abs(t-tilt)));

% rho       = 1.06;   % Density of blood in g/cm^3
% g         = 982;    % Gravitational acceleration in cm/s^2
%                     % Smithsonian physical tables, Smithsonian Institution,
%                     % Seventh Revised Edition, prepared by Frerick E Fowle,
%                     % Vol 71, #4, 1921
% 
% 
% h  =  20; 
% theta     = 60;             % Maximal tilt-angle, in degrees.
% a         = 60/tilt_length;           % Tilt-speed
% conv      = 1333.22;        % conversion cgs (cm,grams,sec) units to mmHg
% arg = zeros(size(t));  
% for z = 1:length(t)
%     if t(z) < tilt
%        arg(z) = 0;
%     elseif t(z) < tilt+tilt_length
%        arg(z) = a.*(t(z)-tilt);
%     else
%         arg(z) = theta; %No tilting down
%     end
% end

% Calculate the hydrostatic pressure in mmHg, 1Pa=0.0075mmHg
% rhogh = rho.*g.*h.*(sin(arg.*pi./180)./conv);

% new_pressure = (p-rhogh-80)*.65+75;%Line up the pulse pressures
new_pressure = (p-80)*.65+75;%Line up the pulse pressures
new_pressure(1:ind_tilt) =new_pressure(1:ind_tilt)+[1:ind_tilt]'.*(-20/ind_tilt)+20;%Take out non-stationary stuff - machine drift

hr = dat(:,3);
new_hr1 = (hr./60-.46); %Shift HR down
new_hr1(ind_tilt:end) = new_hr1(ind_tilt:end)*1.05; %Experience tachycardia faster so it fits in the image
new_hr = new_hr1;
new_hr(ind_tilt+10*Fs:end) = new_hr1(ind_tilt+10*Fs:end).*1.1; %Experience tachycardia faster so it fits in the image

figure(1)
clf
plot(t,new_hr)
xlim([1840-75,1840+75])

writematrix([t,new_hr,new_pressure],'Adjusted_POTS_pp_hr.txt')