% ***********************************************************************************
%                          L O A D  G L O B A L  for
%            C A R D I O V A S C U L A R   S Y S T E M S   M O D E L
% ***********************************************************************************
% Copyright (C) 2019 A.L. Colunga(1), N.P. Woodall(1), M.S. Olufsen(1), and B.E. Carlson
% 
% 1) Department of Mathematics, NC State University 
% 2) Molecular and Integrative Physiology, University of Michigan
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, and merge the Software subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% Any work resulting from the use of the Software should cite the associated 
% article arXiv:1812.11857
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, OR DATA, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
% ***********************************************************************************

% ***********************************************************************************
% This file conatins the right heart catheterization data for Patient 233
% taken ob 29 October 2015.
%
% Output: 
% Patient233_Values  structure with data value
% Patient233_Fields  structure with field names
% data               data structure 
% ***********************************************************************************
function data = Patient
% Patient 233 with right heart catheterization on 29 October 2015
% P_RVsyst  = 40;                             % Systolic RV pressure (mmHg)     
% P_RVdiast = 4;                              % Diastolic RV pressure (mmHg)    
% P_PAsyst  = 21;                             % Systolic pulm art press (mmHg)  
% P_PAdiast = 8;                             % Diast pulm art press (mmHg)     
% P_PV      = 5;                        % Average pulm wedge press (mmHg) 
P_SAsys  = 120;                            % Systolic arterial press (mmHg)    
P_SAdia = 80;                             % Diastolic arterial press (mmHg)    
CO        = 5*1000/60;                      % Cardiac outpt (mL/sec)  
BW        = 75;                           % Body weight (kg)                
Hgt       = 180;                            % Height (cm)                     
Gender    = 2;                              % Gender (1=female, 2=male)       
Ave_HR    = 60;                          % Average heart rate (beats/min)    
% NC        =30;   
% Structure with data value
Patient_Values = {P_SAsys P_SAdia ...
    Ave_HR CO BW Hgt Gender};

% Structure with field names
Patient_Fields = {'P_SAsys' 'P_SAdia' ...
    'Ave_HR' 'CO' 'BW' 'Hgt' 'Gender'};
% Data structure  
data = cell2struct(Patient_Values, Patient_Fields, 2);
clear Patient_Values Patient_Fields
end
