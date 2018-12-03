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

num_simulations = 300;
t_max = 0.01;
%% Define initial conditions
L = 10000;
a0 = 2000;
b0 = 1000;
g0 = 1000;
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


figure
for iter = 1:num_simulations
    iter
    hold on
    memory{iter} = run_simulation(initial_cond, SA, k_array, t_max);
end


%% Define disease initial conditions
L = 100;
a0 = 2000;
b0 = 1000;
g0 = 1000;
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


figure
for iter = 1:num_simulations
    iter
    hold on
    memory_disease{iter} = run_simulation(initial_cond, SA, k_array, t_max);
end



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




%%
t_space = linspace(0,t_max,200);

% hist3(asdf,'CdataMode','auto','nbins',[200,20], 'EdgeColor','none','FaceColor','interp')

[L_WT, L_KO] = paired_time_hist3(memory, memory_disease, 2, t_space);
%%
figure 
subplot(1,2,1)
hist3(L_WT,'CdataMode','auto','nbins',[200,20], 'EdgeColor','none','FaceColor','interp')
subplot(1,2,2)
hist3(L_KO,'CdataMode','auto','nbins',[200,20], 'EdgeColor','none','FaceColor','interp')
%%
hold on
hist3(wam,'CdataMode','auto','nbins',[200,20], 'EdgeColor','none','FaceColor','interp')
hist3(wam2,'CdataMode','auto','nbins',[200,20], 'EdgeColor','none','FaceColor','interp')
function dataset = gimme_the_histogram(memory, var_indx, t_space)
    % Gets histogram data for one dataset.
    times = [];
    datas = [];
    for tp = 1:length(t_space)
        time = t_space(tp);
        temp = get_timept_data(memory, var_indx, time, 1);
        times = [times; temp(:,1)];
        datas = [datas; temp(:,2)];
    end
    dataset = [datas,times];
end


function [dataset1, dataset2] = paired_time_hist3(memory1, memory2, var_indx, t_space)
    % helper function for creating paired disease state and normal state
    % plots using hist3
    t1 = [];
    d1 = [];
    t2 = [];
    d2 = [];
    for tp = 1:length(t_space)
        time = t_space(tp);
        data1 = get_timept_data(memory1, var_indx, time, 1);
        data2 = get_timept_data(memory2, var_indx, time, 1);
        t1 = [t1; data1(:,1)];
        t2 = [t2; data2(:,1)];
        
        d1 = [d1; data1(:,2)];
        d2 = [d2; data2(:,2)];
    end
    dataset1 = [d1,t1];
    dataset2 = [d2,t2];
end
        

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
