function main
%% Main body of simulation
    clear; close all; clc; tic;
    addpath(genpath('functions'))
    % Initialize Parameters
    e = [10,1000,10];
    s = [100,100,5];
    es = [0,0,0];
    p = [0,0,0];
    c = [0.1; 1; 1];
    SA = [1 0 0 ; 1 -1 1 ; 0 -1 1 ; -1 1 -1];
    t_max = [40, 10, 10];

    simulations = 50;
    memory = cell(simulations,length(t_max));
    memory_approx = cell(simulations,length(t_max));

    % Run <iter> number of simulations for each case <run> 
    for iter = 1:simulations
        exception = zeros(1,2);
        for run = 1:length(t_max)
            memory{iter,run} = run_simulation(p(run),e(run),s(run),es(run),c,SA,t_max(run),iter,run);
            memory_approx{iter,run} = run_approx_sim(p(run),e(run),s(run),es(run),c,SA,t_max(run),iter,run);
        end
    end
    fprintf('Gillespie Simulation Completed\t');
    toc
    % Removing extraneous data (Z, alpha, rm)
    cleaned_data = rm_extras(memory,2,10);
    cleaned_approx_data = rm_extras(memory_approx,3,5);
    % Recovering only End state data
    end_states = end_state(cleaned_data);
    end_approx_states = end_state(cleaned_approx_data);
    end_P = collect_var_results(2,end_states);
    end_P_approx = collect_var_results(2,end_approx_states);
    fprintf('Data Processing Completed\t');
    toc
    %% Plot End state distribution
    for p = 1: length(e)
        subplot(length(e),1,p)
        hold on;
        histogram(end_P{p},'Normalization','pdf')
        histogram(end_P_approx{p},'Normalization','pdf')
        title(sprintf('e = %i, s = %i, # sim = %i',e(p),s(p),simulations));
        xlabel('units of P');
        ylabel('proportion');
        legend('Exact','Approximation');
        fprintf('(e,s,n_sim,t_max) = (%i,%i,%i,%i)\n',e(p),s(p),simulations,t_max(p));
        fprintf('\tExact PMF Mean = %f\n',mean(end_P{p}));
        fprintf('\tExact PMF Mode = %f\n',mode(end_P{p}));
        fprintf('\tExact PMF Median = %f\n',median(end_P{p}));
        fprintf('\tApproximate PMF Mean = %f\n',mean(end_P_approx{p}));
        fprintf('\tApproximate PMF Mode = %f\n',mode(end_P_approx{p}));
        fprintf('\tApproximate PMF Median = %f\n',median(end_P_approx{p}));
    end
    
    
    
    filename = 'workspace_data.mat';
    save(filename);
end
    

