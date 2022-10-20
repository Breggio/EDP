%% Estimation of a site where motion of a moving vehicle started using radar data
% Coordinate transformation of measurements 

% Written by Irina Yareshko and Luca Breggion, Skoltech 2022

close all 
clear
clc

set(0,'defaulttextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%

A = importdata('z1.txt');