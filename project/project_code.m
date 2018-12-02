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

k_aL            = 7.8e-2;
k_r_aL          = 0.4;
k_bL            = 7.9e-2;
k_r_bL          = 0.22;

k_aLb           = 1.5e-2;               % 5e-6:1.5e-2
k_r_aLb         = 1.7e-3;             % 1.7e-3:0.22
k_bLg           = 1.5e-2;               % 5e-6:1.5e-2
k_r_bLg         = 1.7e-3;             % 1.7e-3:4
k_bLa           = 1.5e-2;               % 5e-6:1.5e-2
k_r_bLa         = 1.7e-3;             % 1.7e-3:0.6
k_abLg          = 1.5e-2;               % 5e-6:1.5e-2
k_r_abLg        = 1.7e-3;             % 1.7e-3:4
k_bgLa          = 1.5e-2;               % 5e-6:1.5e-2
k_r_bgLa        = 1.7e-3;             % 1.7e-3:0.6

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

%% Define initial conditions
L = 1000;
a0 = 500;
b0 = 500;
g0 = 500;
% Set IC
initial_cond = zeros(num_reagents, 1);
initial_cond(1) = L;
initial_cond(2) = a0;
initial_cond(3) = b0;
initial_cond(4) = g0;

initial_cond(19) = initial_cond(2);
initial_cond(20) = initial_cond(3);
initial_cond(21) = initial_cond(4);

%%
hold on
for i = 1:1
    i
    Z = zeros(num_rxn, 1);
    alpha = zeros(num_rxn, 1);
    t = 0;
    identity = eye(num_rxn);

    num_simulations = 5;
    t_max = 2.5;

    % history = [t, alpha', z', zeros(1,num_rxn), initial_cond'];
    x = initial_cond;

    %%%%%% start of loops
    rxn_order= [];
    history = [t,alpha',Z',zeros(1,num_rxn),x'];
    while t<t_max
        alpha = update_alpha(Z, SA, k_array, initial_cond);
        rm = alpha/sum(alpha);
        rxn_q = rxn_number(rm);
        rxn_order = [rxn_order; rxn_q];
        try
            Z = Z + identity(:,rxn_q);
        catch
    %         fprintf('Rxn terminated early at run %i, iter %i\n'
            fprintf('RXN terminated early');
            break;
        end
        alpha;
        tau = exprnd(1/sum(alpha));
        if isnan(tau)
            tau = 0;
        end

        t = t + tau;
    %     w = waitforbuttonpress
        x_advance = initial_cond + SA*Z;
        history = [history;t,alpha',Z',rm',x_advance'];
    end
    plot(history(:,1),history(:,77:94)) 
%     plot(history(:,1),history(:,85)) 
end

%%

test = run_simulation(initial_cond, SA, k_array, 2.5);
plot(test(:,1),test(:,77:94)) 
function history = run_simulation(initial_conditions, S_array, K_array, t_max)
    [num_reagents, num_rxn] = size(S_array);
    Z = zeros(num_rxn, 1);
    alpha = zeros(num_rxn, 1);
    identity = eye(num_rxn);
    t = 0;
    x0 = initial_conditions;
    rxn_order= [];
    history = [t,alpha',Z',zeros(1,num_rxn),x0'];
    while t<t_max
        alpha = update_alpha(Z, S_array, K_array, x0);
        rm = alpha/sum(alpha);
        rxn_q = rxn_number(rm);
        rxn_order = [rxn_order; rxn_q];
        try
            Z = Z + identity(:,rxn_q);
        catch
            fprintf('RXN terminated early');
            break;
        end
        tau = exprnd(1/sum(alpha));
        if isnan(tau)
            tau = 0;
        end
        t = t + tau;
        x_advance = x0 + S_array*Z;
        history = [history;t,alpha',Z',rm',x_advance'];
    end
end

function alpha = update_alpha(Z_array, S_array, K_array, initial_conditions)
    [num_reagents, num_rxn] = size(S_array);

    X_updated = S_array * Z_array + initial_conditions;
    reactants_matrix = double(S_array == -1);
    
    for col=1:num_rxn
        reactants_matrix(:,col) = reactants_matrix(:,col) .* X_updated;
    end
    temp = reactants_matrix;
    temp(~reactants_matrix) = 1;
    h_array = prod(temp)';
    alpha = h_array .* K_array;
    
end
    

function quotient = rxn_number(rm)
    r1 = rand;
    sum_R = 0;
    quotient = 0;
    for reaction_quotient = 1:length(rm)
        sum_R = sum_R + rm(reaction_quotient);
        if r1 < sum_R
            quotient = reaction_quotient;
            break
        end
    end
    return
end