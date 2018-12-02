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

num_simulations = 20;

%% Define initial conditions
L = 10000;
a0 = 2000;
b0 = 1000;
g0 = 800;
% Set IC
initial_cond = zeros(num_reagents, 1);
initial_cond(1) = L;
initial_cond(2) = a0;
initial_cond(3) = b0;
initial_cond(4) = g0;

initial_cond(19) = initial_cond(2);
initial_cond(20) = initial_cond(3);
initial_cond(21) = initial_cond(4);


memory = cell(num_simulations,1);
t_max = 0.008;

figure
for iter = 1:num_simulations
    iter
    hold on
    memory{iter} = run_simulation(initial_cond, SA, k_array, t_max);
end


%% Define disease initial conditions
L = 100;
a0 = 2000;
b0 = 500;
g0 = 600;
% Set IC
initial_cond = zeros(num_reagents, 1);
initial_cond(1) = L;
initial_cond(2) = a0;
initial_cond(3) = b0;
initial_cond(4) = g0;

initial_cond(19) = initial_cond(2);
initial_cond(20) = initial_cond(3);
initial_cond(21) = initial_cond(4);


memory_disease = cell(num_simulations,1);
t_max = 0.1;

figure
for iter = 1:num_simulations
    iter
    hold on
    memory_disease{iter} = run_simulation(initial_cond, SA, k_array, t_max);
end


%%

%% REally ugly plotting code 
% Retrieve 
L_data = get_var_time_data(memory,1,2);
a_data = get_var_time_data(memory,1,3);
b_data = get_var_time_data(memory,1,4);
g_data = get_var_time_data(memory,1,5);
aL_data = get_var_time_data(memory,1,6);
bL_data = get_var_time_data(memory,1,7);
abL_data = get_var_time_data(memory,1,8);
bgL_data = get_var_time_data(memory,1,9);
abgL_data = get_var_time_data(memory,1,10);
L_i_data = get_var_time_data(memory,1,11);
a_i_data = get_var_time_data(memory,1,12);
b_i_data = get_var_time_data(memory,1,13);
g_i_data = get_var_time_data(memory,1,14);
aL_i_data = get_var_time_data(memory,1,15);
bL_i_data = get_var_time_data(memory,1,16);
abL_i_data = get_var_time_data(memory,1,17);
bgL_i_data = get_var_time_data(memory,1,18);
abgL_i_data = get_var_time_data(memory,1,19);

L_data_ES = cell2mat(end_state(L_data));
a_data_ES = cell2mat(end_state(a_data));
b_data_ES = cell2mat(end_state(b_data));
g_data_ES = cell2mat(end_state(g_data));
aL_data_ES = cell2mat(end_state(aL_data));
bL_data_ES = cell2mat(end_state(bL_data));
abL_data_ES = cell2mat(end_state(abL_data));
bgL_data_ES = cell2mat(end_state(bgL_data));
abgL_data_ES = cell2mat(end_state(abgL_data));
L_i_data_ES = cell2mat(end_state(L_i_data));
a_i_data_ES = cell2mat(end_state(a_i_data));
b_i_data_ES = cell2mat(end_state(b_i_data));
g_i_data_ES = cell2mat(end_state(g_i_data));
aL_i_data_ES = cell2mat(end_state(aL_i_data));
bL_i_data_ES = cell2mat(end_state(bL_i_data));
abL_i_data_ES = cell2mat(end_state(abL_i_data));
bgL_i_data_ES = cell2mat(end_state(bgL_i_data));
abgL_i_data_ES = cell2mat(end_state(abgL_i_data));


%% Plot block
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(L_data_ES(:,2), 'Normalization','pdf');
% histogram(L_i_data_ES(:,2), 'Normalization','pdf');
title("Ligand [L] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "L_ES.png")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(a_data_ES(:,2), 'Normalization','pdf');
% histogram(a_i_data_ES(:,2), 'Normalization','pdf');
title("Alpha subunit [a] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "a_ES.png")


figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(b_data_ES(:,2), 'Normalization','pdf');
% histogram(b_i_data_ES(:,2), 'Normalization','pdf');
title("Beta subunit [b] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "b_ES.png")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(g_data_ES(:,2), 'Normalization','pdf');
% histogram(g_i_data_ES(:,2), 'Normalization','pdf');
title("Gamma subunit [g] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "g_ES.png")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(aL_data_ES(:,2), 'Normalization','pdf');
% histogram(aL_i_data_ES(:,2), 'Normalization','pdf');
title("Alpha-Ligand subunit [aL] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "aL_ES.png")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(bL_data_ES(:,2), 'Normalization','pdf');
% histogram(bL_i_data_ES(:,2), 'Normalization','pdf');
title("Beta-Ligand subunit [bL] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "bL_ES.png")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(abL_data_ES(:,2), 'Normalization','pdf');
% histogram(abL_i_data_ES(:,2), 'Normalization','pdf');
title("Alpha-Beta-Ligand subunit [abL] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "abL_ES.png")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(bgL_data_ES(:,2), 'Normalization','pdf');
% histogram(bgL_i_data_ES(:,2), 'Normalization','pdf');
title("Beta-Gamma-Ligand subunit [bgL] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "bgL_ES.png")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
histogram(abgL_data_ES(:,2), 'Normalization','pdf');
% histogram(abgL_i_data_ES(:,2), 'Normalization','pdf');
title("Alpha-Beta-Gamma-Ligand subunit [bgL] Count at 5s");
xlabel('Units [count]');
ylabel('Probability');
% legend('External','Internal');
saveas(gcf, "abgL_ES.png")



%%
get_timept_data(memory, 2, 0.008, 1)

function timept_data = get_timept_data(memory_matrix, variable_indx, time_pt, time_indx)
    % Given a history matrix and a variable ID number, gets the value for
    % that variable at the time_pt across all simulations. time_indx
    % indicates location of the time array.
    variable_data = get_var_time_data(memory_matrix, time_indx, variable_indx);
    if time_pt == inf
        timept_data = cell2mat(end_state(variable_data));
    else
        [num_sims, num_runs] = size(variable_data);
        timept_data = zeros(num_sims, 2);
        for i = 1:num_sims
            time_list = variable_data{i}(:,1);
            for t = 1:length(time_list)
                time = variable_data{i}(t,1);
                data = variable_data{i}(t,2);
                if time > time_pt
                    timept_data(i,1) = time;
                    timept_data(i,2) = data;
                    break
                else
                    continue
                end
            end
        end
    end
end
    

function var_data = get_var_time_data(memory_matrix, time_indx, variable_indx)
    % Gets all of the time data for a given variable.
    [num_sims, num_runs] = size(memory_matrix);
    var_data = cell(num_sims, num_runs);
    for i = 1:num_sims
        extraction = memory_matrix{i}(:,variable_indx);
        time_extract = memory_matrix{i}(:,time_indx);
        merge = [time_extract, extraction];
        var_data{i,num_runs} = merge;
    end
end


function history = run_simulation(initial_conditions, S_array, K_array, t_max)
    % Runs gillespie algorithm given a set of initial conditions,
    % stochiometry matrix, reaction rate coefficients, and maximum time to
    % run the simulation.
    [num_reagents, num_rxn] = size(S_array);
    Z = zeros(num_rxn, 1);
    alpha = zeros(num_rxn, 1);
    identity = eye(num_rxn);
    t = 0;
    x0 = initial_conditions;
    rxn_order= [];
%     history = [t,alpha',Z',zeros(1,num_rxn),x0'];
    history = [t,x0'];
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
            break
        end
        t = t + tau;
        x_advance = x0 + S_array*Z;
%         history = [history;t,alpha',Z',rm',x_advance'];
        history = [history;t,x_advance'];
    end
    plot(history(:,1),history(:,2:8));
end


function alpha = update_alpha(Z_array, S_array, K_array, initial_conditions)
    % Updates the propensity matrix and determines the reaction that will
    % occur based on the number of elements per reactant and the
    % probability of each reaction occuring.
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
    % Randomly determines the reaction that will occur by given the
    % probability distribution of a reaction occuring.
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


function clean_data = rm_extras(memory,start,stop)
    % function removes redundant columns from dataset, takes in start and
    % stop range for removing columns by index.
    [iter, run] = size(memory);
    clean_data = cell(iter,run);
    for i = 1:iter
        for r = 1:run
            extract = memory{i,r};
            extract(:,start:stop)=[];
            clean_data{i,r} = extract;
        end
    end
    return
end


function end_results = end_state(cleaned_data)
    % Function grabs the last element in a simulation, takes the last
    % possible occurence in the history of the simulation (last time point)
    [iter, run] = size(cleaned_data);
    end_results = cell(iter,run);
    for i = 1:iter
        for r = 1:run
            extract = cleaned_data{i,r}(end,:);
            end_results{i,r} = extract;
        end
    end
    return
end
