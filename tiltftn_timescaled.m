% This function returns the hydrostatic pressure, calculated based upon the
% time elapsed in the simulation.

function [rhogh,arg] = tiltftn_timescaled(t,pars)  % ADD TD (Tilt Down) to function


% Define tilt function parameters
rho       = 1.06;   % Density of blood in g/cm^3
g         = 982;    % Gravitational acceleration in cm/s^2
                    % Smithsonian physical tables, Smithsonian Institution,
                    % Seventh Revised Edition, prepared by Frerick E Fowle,
                    % Vol 71, #4, 1921

ts =  pars(end-2);
td =  pars(end-1);
h  =  pars(end);

tdur      = td - ts;        %tilt end
theta     = 70;             % Maximal tilt-angle, in degrees.
a         = 60/14;           % Tilt-speed
conv      = 1333.22;        % conversion cgs (cm,grams,sec) units to mmHg
  
if t < ts;
   arg = 0;
elseif t < ts+14
   arg = a.*(t-ts);
elseif t < td;
    arg = theta;
elseif t < td+14
    arg = theta - a*(t-td);
else
    arg = 0;
end;


% Calculate the hydrostatic pressure in mmHg, 1Pa=0.0075mmHg
rhogh = rho*g*h*(sin(arg*pi/180)/conv); 
end

