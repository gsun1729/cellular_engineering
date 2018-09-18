%This assignment is Cell Engineering HW#1 for Fall 2018
%Author:  Gordon Sun
%Date:    20180915
%
%This program uses nonlinear least squares on wild-type and mutant protein/ligand interaction data to...
%Assumptions:  Assuming steady state.
%This file outputs 4 graphs

function myFlag = CellEng_HW1_2018Template()
clear; close all; clc;
%Wild-Type Data.  [L] (nM) vs signal

L0   = [1 2 4 8 10 20 40 80 100 200 400 800 1600 3200 6400 1e4];   % Given ligand concentration (independent variable)
WT   = [0.0178
    0.0441
    0.0826
    0.2060
    0.2222
    0.336
    0.5533
    0.6457
    0.7207
    0.8811
    0.9396
    0.9181
    0.9686
    0.9892
    0.9846
    0.9965
    ]';   % WT Signal (dependent variable)
MUT  = [0.0312
    0.0786
    0.1454
    0.2432
    0.2727
    0.3208
    0.4602
    0.4346
    0.6288
    0.5175
    0.5519
    0.5978
    0.5546
    0.7384
    0.7198
    0.6828
    ]';   % Mut Signal (dependent variable)

beta0 = [0, 0];      % Parameter 1:  Initial Guesses: Constant of Proportionality between fraction saturated and signal
% Parameter 2:  Initial Guesses: Kd [=] nM

[WT_beta, WT_res] = nlinfit(L0, WT, @binding, beta0);    % WT_beta = [proportionality constant, Kd] , ** This needs to be filled in to run **

%Plots experimental wild-type data
figure;
semilogx(L0, WT);  % what you are plotting
xlabel('Ligand Concentration (nM)', 'FontSize', 16);
ylabel('WT Signal', 'FontSize', 16);
title('WT Protein Signal vs. Ligand Concentration', 'FontSize', 18);
axis([ ]);

%Plots fit of wild-type data
hold on;
p =  L0; % this is your array of ligand points to use for fitting
WT_fit = WT_beta(1) * p ./ (p + WT_beta(2)); % this will be a function of WT_beta(1) WT_beta(2), and p
semilogx(L0, WT_fit); % what you are plotting
legend('Raw WT data', 'WT Fit');
hold off;

%Plots Residuals
figure;
semilogx(L0, WT_res);  % what you are plotting
xlabel('Ligand Concentration (nM)', 'FontSize', 16);
ylabel('Residuals', 'FontSize', 16);
title('WT Residuals', 'FontSize', 18);
axis([  ]);

%Repeat with Mutant Data Below.  [L] (nM) vs Mutant Signal

return; 


%------------------Function that describes the curve that data points are fit to--------------
function signal = binding(beta, L)
Var1   = beta(1);
Var2   = beta(2);
signal  =  Var1 * L ./ (L + Var2); % Will be a function of L and what you name as Variables (Var1 and Var2), equation from class
return;