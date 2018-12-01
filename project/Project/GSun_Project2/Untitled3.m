clear; close all; clc; tic;

%% Initial parameters initialization
SA = [
    -1	1	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    -1	1	0	0	0	0	0	0	-1	1	0	0	-1	1	0	0	0	0	0	-1	0	0	1	0	0
    0	0	-1	1	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	1	0
    0	0	0	0	0	0	-1	1	0	0	-1	1	0	0	0	0	0	0	0	0	0	-1	0	0	1
    1	-1	0	0	-1	1	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0
    0	0	1	-1	0	0	-1	1	-1	1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0
    0	0	0	0	1	-1	0	0	1	-1	-1	1	0	0	0	0	-1	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	1	-1	0	0	0	0	-1	1	0	0	0	-1	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	1	-1	1	-1	0	0	0	0	-1	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
];

[num_reagents, num_rxn] = size(SA);

k_aL            = 7.8e6;
k_r_aL          = 0.4;
k_bL            = 7.9e5;
k_r_bL          = 0.22;

k_aLb           = 1.5e-2;               % 5e-6:1.5e-2
k_r_aLb         = 0.22;             % 1.7e-3:0.22
k_bLg           = 1.5e-2;               % 5e-6:1.5e-2
k_r_bLg         = 4;             % 1.7e-3:4
k_bLa           = 1.5e-2;               % 5e-6:1.5e-2
k_r_bLa         = 0.6;             % 1.7e-3:0.6
k_abLg          = 1.5e-2;               % 5e-6:1.5e-2
k_r_abLg        = 4;             % 1.7e-3:4
k_bgLa          = 1.5e-2;               % 5e-6:1.5e-2
k_r_bgLa        = 0.6;             % 1.7e-3:0.6

k_int           = (7e-3) / 60;
k_sig_int       = (4e-2) / 60;

k_array = [
    k_aL
    k_r_aL
    k_bL
    k_r_bL
    k_aLb
    k_r_aLb
    k_bLg
    k_r_bLg
    k_bLa
    k_r_bLa
    k_abLg
    k_r_abLg
    k_bgLa
    k_r_bgLa
    k_int
    k_int
    k_int
    k_sig_int
    k_sig_int
    k_int
    k_int
    k_int
    k_int
    k_int
    k_int
    ];

%% Define intial conditions
L = 100;
a0 = 100;
b0 = 100;
g0 = 100;
%% Set IC
initial_cond = zeros(num_reagents, 1);
initial_cond(1) = L;
initial_cond(2) = a0;
initial_cond(3) = b0;
initial_cond(4) = g0;

initial_cond(19) = initial_cond(2);
initial_cond(20) = initial_cond(3);
initial_cond(21) = initial_cond(4);

%%
Z = zeros(num_rxn, 1);
alpha = zeros(num_rxn, 1);
t = 0;
identity = eye(num_rxn);

num_simulations = 5;
t_max = 50;

history = [t, alpha', z', zeros(1,num_rxn), initial_cond'];


%%%%%% start of loops

S_Z = SA * Z;

    

function alpha = update_alpha(initial_cond, k_array, S_Z)
    alpha = zeros(length(k_array),1);
    
end

% 
% function history = run_simulation(initial_cond, stochi_matrix, rxn_coeff, t_max)
%     
% end