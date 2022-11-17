clc; close all; clear

%% Read Data

pathRoot = 'C:\OneDrive - Delft University of Technology\3. Semester - Studienunterlagen\Engineering Geophysics\assignment 2 Downhole seismic';

cptData = readtable(fullfile(pathRoot,'downhole-cptRead2Matlab.txt'));

tvtData = readtable(fullfile(pathRoot,'downhole-traveltimes.xls'));

%% DHT Calculate true distance and interval velocity
% Define Parameters

Eg = 0.6;    % elevation of the top of the geophone hole [m]
Es = 0;      % elevation of the ground surface in contact with the energy source at the center of the energy source
X  = 1.5;    % horizontal distance between the center of the energy source and the geophone borehole/sounding


% The following equation determines the straight-line
% slant distance, LR, from source to geophone 
tvtData.Lr = sqrt((Es-Eg+tvtData.z_m_).^2 + X^2);

intervalDataDHT    = table;
intervalDataDHT.dR = tvtData.Lr(2:end) - tvtData.Lr(1:end-1);
intervalDataDHT.dT = (tvtData.t_s_(2:end) - tvtData.t_s_(1:end-1));

intervalDataDHT.Vs =  intervalDataDHT.dR ./ intervalDataDHT.dT;

% prepare for plot
plotDataDht    = table;
plotDataDht.dZ = [tvtData.z_m_(1)-Eg; diff(tvtData.z_m_)];
plotDataDht.Z  = cumsum(plotDataDht.dZ);
plotDataDht.Vs = [intervalDataDHT.Vs;NaN];
%% SCPT Calculate true distance and interval velocity
% Cone resistance, qc The force acting on the cone, Qc, divided by the projected area of the cone, Ac. 
%  qc = Qc / Ac

% Sleeve friction, fsThe frictional force acting on the friction sleeve, Fs, divided by its surface area, As. 
%  fs = Fs / As
% 
% "Depth";"qc; "fs"; „Inclination"; „dummy";  "fs/qc"; f
% "[m]";"[MN/m²]";"[MN/m²]";"[°]", "[]";"[]"


% Empirical relatioship Hegazy and Mayne (1995) - Sand
cptData.Vs = 12.02.*(cptData.qc*1000).^(0.319).*(cptData.fs*1000).^(-0.0466);

% Empirical relatioship Hegazy and Mayne (1995) - Sand
cptData.Vs1 = 1.75*(cptData.qc*1000).^(0.627);

% prepare for plot
% plotDataCpt    = table;
% plotDataCpt.dZ = [0; diff(cptData.Tiefe)];
% plotDataCpt.Z  = cumsum(plotDataCpt.dZ);
% plotDataCpt.Vs = [intervalDataCPT.Vs;NaN];
%% Plot Results


[zPlotDHT,vsPlotDHT] = stairs(plotDataDht.Z,plotDataDht.Vs);
[zPlotCPT,vsPlotCPT] = stairs(cptData.Depth, cptData.Vs);
[zPlotCPT1,vsPlotCPT1] = stairs(cptData.Depth, cptData.Vs1);


hold on
plot(vsPlotDHT,zPlotDHT)
plot(vsPlotCPT,zPlotCPT)
plot(vsPlotCPT1,zPlotCPT1)

set(gca, 'YDir','reverse')
ytickformat('%g m')
xtickformat('%g m/s')
title('Downhole Seismic Test (DHT)')
xlabel('interval velocity')
ylabel('depth')