%% The relationship between solar radio flux F10.7 and sunspot number 

close all; clearvars; clc;
set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1) Download monthly mean sunspot number and solar radio flux F10.7cm

% Format of data: 
%1 column - year  
%2 column – month 
%3 column – monthly solar radio flux at 10.7 cm 
%4 column – monthly sunspot number

data = load('data_group8.mat');


%% 2) Plot the monthly mean sunspot number and solar radio flux F10.7 cm 
% for visual representation. 




