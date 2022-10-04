%% The relationship between solar radio flux F10.7 and sunspot number 

close all; clear all; clc;

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% 1) Download monthly mean sunspot number and solar radio flux F10.7cm

data = load('data_group8.mat');
y = data.data_group8(:,1); % year
m = data.data_group8(:,2); % month
m_flux =  data.data_group8(:,3); % monthly solar radio flux at 10.7 cm 
m_sun =  data.data_group8(:,4); % monthly sunspot number

%% 2) Plot the monthly mean sunspot number and solar radio flux F10.7 cm 
% for visual representation. 

%% 3) Make scatter plot between monthly mean sunspot number and solar radio flux F10.7 cm

figure(1)


