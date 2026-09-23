close all
clear all

load("processed_data_classification_script_50msbin.mat","normall","LA","MLA",'highNLA')

XLA=normall(:,[LA;MLA]); % normalized data to mechanical trials only: Mechanical LinearlyResponsive neurons
XNN=normall(:,highNLA); % normalized data to mechanical trials only: Mechanical Non NonLinearly Responsive neurons

% Paramètres
Fs_original = 20;   % original Frequency (Hz)
Fs_target = 200;     % targeted Frequency (Hz)
pressure = [1:20];
pressure_target = 1:0.1:20; 

y_lisse = interpUpsampleSmooth(pressure, XLA, pressure_target, 'spline', 20);
y_lisse1 = interpUpsampleSmooth(pressure, XNN, pressure_target, 'spline', 20);

figure; 
subplot(2,2,1)  % raw data LA
plot(XLA)
title('raw data LA')
xlabel('force')
ylabel('perc. of activation')
subplot(2,2,2)  % smoothed data LA
plot(y_lisse)
title('smoothed data LA')
xlabel('force')
ylabel('perc. of activation')
subplot(2,2,3)  % raw data NLA
plot(XNN)
title('raw data NLA')
xlabel('force')
ylabel('perc. of activation')
subplot(2,2,4)  % smoothed data NLA
plot(y_lisse1)
title('smoothed data NLA')
xlabel('force')
ylabel('perc. of activation')



