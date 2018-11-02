

function y0 = HW4_Template2018()
clear all;  close all; clc;

%---------------------------------------------------SET UP--------------------------------------------------------------

%Set Parameters
B       = 1e-7;             % 1 / person day
c       = 50;               % 
p       = 0.5;              % 
N       = 1e7;             % 
r1      = 1/3;              % 1 / days
r2      = 1/8;              % 1 / days
r3      = 1/3;              % 1 / days
r4      = 1/12;             % 1 / days
n       = 5000;             % 
u       = 50;               % 1 / day (TV)
h       = 0.9;              % 
alpha   = 1/5;            % 1 / days
v0      = 0.975;            % 
v1      = 0.975;           % 
del     = 0.3;             % 
f       = 1e-6;             % 
I_t0    = 1e3;
 
%Pack Parameters
theta = [B, c, p, N, r1, r2, r3, r4, n, u, h, alpha, v0, v1, del, f]; 
% Set initial conditions for Phase I
tspan = [0 5];
initial_conditions = [N I_t0 0 0 0 0 0];


% Solve ODEs for Phase I
[t, y]  = ode23s(@phaseI, tspan, initial_conditions,  [], theta);
final_state = y(end,:);
pad = zeros(length(y(:,1)), 10);
y_padded = [y(:,1:2),pad, y(:,3:end)];
%------------------------------------------------Trace Vaccination----------------------------------------------------------------

% Set initial conditions for Phase II TV
p2_initial_cond = [final_state(1:end-2) 0 0 0 0 0 0 0 0 0 0 final_state(end-1:end)];
tspan2 = [0 300];
% % Variables
% % S0      = y(1);
% % I0_1    = y(2);
% % I0_2    = y(3);
% % I0_3    = y(4);
% % I0_4    = y(5);
% % Q_0     = y(6);
% % Q_1     = y(7);
% % Q_2     = y(8);
% % Q_3     = y(9);
% % H       = y(10);
% % S1      = y(11);
% % I1_1    = y(12);
% % I1_2    = y(13);
% % I1_3    = y(14);
% % I1_4    = y(15);
% % Z       = y(16);
% % D       = y(17)
% Solve ODEs for Phase II TV

[t2, y2]  = ode23s( @phaseII, tspan2, p2_initial_cond, [], theta);

% Keep track of infected, queued, and dead

t2_inc = t2 + 5;
TV_merged_time = [t; t2_inc(2:end)];
merged_TV = [y_padded; y2(2:end,:)];
% TV_I = sum(merged_TV(:,2:5), 2) + sum(merged_TV(:,12:15), 2);
TV_I = sum(merged_TV(:,2:5), 2) + sum(merged_TV(:,12:15), 2);
TV_Q = sum(merged_TV(:,6:9), 2);
TV_D = merged_TV(:,17);

% %------------------------------------------------Mass Vaccination-------------------------------------------------------------------
% 
% Set initial conditions for Phase II MV
MV_initial_cond = [0 0 0 0 final_state(5) final_state(1) final_state(2) final_state(3) final_state(4) 0 0 0 0 0 0 final_state(end-1:end)];
theta(10) = 200;
tspan3 = [0 50];

%Solve ODEs for Phase II MV

[t3, y3]  = ode23s(@phaseII, tspan3, MV_initial_cond,  [], theta);
% Keep track of infected, queued, and dead
t3_inc = t3 + 5;
MV_merged_time = [t; t3_inc(2:end)];
% plot(t3, y3(:,))
merged_MV = [y_padded; y3(2:end,:)];
MV_I = sum(merged_MV(:,2:5), 2) + sum(merged_MV(:,12:15), 2) + sum(merged_MV(:,7:9), 2);
MV_Q = sum(merged_MV(:,6:9), 2);
MV_D = merged_MV(:,17);

% return
%-------------------------------------------------------------------------Plots------------------------------------------------

figure;
[AX,H1,H2] = plotyy(TV_merged_time, TV_Q, TV_merged_time, TV_I);
set(get(AX(1),'Ylabel'),'String','Number in Queue');
set(get(AX(2),'Ylabel'),'String','Number Infected');
set(H1,'LineStyle',':');
set(H2,'LineStyle','-');
title('Trace Vaccination');
xlabel('Time (days)');
legend(H1, 'queue');
legend(H2, 'infected');


figure;
[AX,H3,H4] = plotyy(MV_merged_time, MV_Q, MV_merged_time, MV_I);
set(get(AX(1),'Ylabel'),'String','Number in Queue');
set(get(AX(2),'Ylabel'),'String','Number Infected');
set(H3,'LineStyle',':');
set(H4,'LineStyle','-');
title('Mass Vaccination');
xlabel('Time (days)');
legend(H3, 'queue');
legend(H4, 'infected');


%------------------------------------------------Death Rate Comparison------------------------------------------------------------    
    
% CDC policy (CDC guidelines are modeled by switching from TV to MV 28 days after the start of TV (e.g. on day 33 in the base case))


% Set initial conditions 
CDC_initial_cond = [final_state(1:end-2) 0 0 0 0 0 0 0 0 0 0 final_state(end-1:end)];
tspan_TV = [0 28];
theta(10) = 50;
% Solve
[t4, y4]  =  ode23s( @phaseII, tspan_TV, CDC_initial_cond, [], theta);

t4_inc = t4 + 5;
CDC_merged_time = [t; t4_inc(2:end)];
CDC_merged = [y_padded; y4(2:end,:)];
% Set initial conditions 
TV_end_cond = CDC_merged(end, :);

CDC_MV_IC = TV_end_cond;
adder = zeros(size(TV_end_cond));
adder(1:4) = TV_end_cond(1:4) * -1;
adder(6:9) = TV_end_cond(6:9);
CDC_MV_IC = CDC_MV_IC + adder;


tspan_MV = [0 100];
theta(10) = 200;
% Solve
[t5, y5]  = ode23s( @phaseII, tspan_MV, CDC_MV_IC, [], theta);

t5_inc = t5 + 33;
CDC_merged_time = [CDC_merged_time; t5_inc(2:end)];
CDC_merged = [CDC_merged; y5(2:end, :)];

% Plot TV vs MV vs CDC Policy

figure;

plot(TV_merged_time,TV_D);
hold on;
plot(MV_merged_time, MV_D, 'r');
plot(CDC_merged_time, CDC_merged(:,end), 'g');
title('Small Pox Deaths');
xlabel('Time (Days)');
ylabel('Deaths (count)');
legend('TV', 'MV', 'CDC Policy');
set(gca, 'YScale', 'log')
hold off;




%------------------------Sensitivity Analysis----------------------------

% Fraction of population initially infected

tspan   = [0 5];
tspan2  = [0 400];
dead    = zeros(50,1);
dead2   = zeros(50,1);
frac    = zeros(50,1);
for f=1:50
    frac(f) = f / 5000; % such as f divided by something 
    i0      = frac(f) * N;% such as frac(f) multiplied by something
    init    = [N i0 0 0 0 0 0];
    [t, y]  = ode23s(@phaseI, tspan, init, [], theta);
    final_state = y(end, :);
    
    init2   = [final_state(1:end-2) 0 0 0 0 0 0 0 0 0 0 final_state(end-1:end)];
    theta(10) = 50;
    [t2, y2] = ode23s(@phaseII, tspan2, init2, [], theta);
    
    dead(f) = y2(end, 17);    
    
    init4   = [0 0 0 0 final_state(5) final_state(1) final_state(2) final_state(3) final_state(4) 0 0 0 0 0 0 final_state(end-1:end)];
    theta(10) = 200;
    [t4, y4]  = ode23s(@phaseII, tspan2, init4, [], theta);
   
    dead2(f)  = y4(end, 17);
end

figure;
plot(frac, dead);
hold on;
plot(frac, dead2, 'r');
title(' ');
xlabel('Fraction of Population Initially Infected');
ylabel('Number of Deaths');
legend('TV', 'MV');
hold off;



% Choose another parameter to vary to conduct sensitivity analysis
% Fraction of vacinnators in the population

tspan   = [0 5];
tspan2  = [0 400];
dead    = zeros(50,1);
dead2   = zeros(50,1);
frac    = zeros(50,1);
i0 = 1000;
for f=1:50
    theta(9) = f * 100;
    frac(f) = f / N; % such as f divided by something 
    init    = [N i0 0 0 0 0 0];
    [t, y]  = ode23s(@phaseI, tspan, init, [], theta);
    final_state = y(end, :);
    
    init2   = [final_state(1:end-2) 0 0 0 0 0 0 0 0 0 0 final_state(end-1:end)];
    theta(10) = 50;
    [t2, y2] = ode23s(@phaseII, tspan2, init2, [], theta);
    
    dead(f) = y2(end, 17);    
    
    init4   = [0 0 0 0 final_state(5) final_state(1) final_state(2) final_state(3) final_state(4) 0 0 0 0 0 0 final_state(end-1:end)];
    theta(10) = 200;
    [t4, y4]  = ode23s(@phaseII, tspan2, init4, [], theta);
   
    dead2(f)  = y4(end, 17);
end

figure;
plot(frac, dead);
hold on;
plot(frac, dead2, 'r');
title(' ');
xlabel('Number of Vaccinators per capita');
ylabel('Number of Deaths');
legend('TV', 'MV');
set(gca, 'YScale', 'log')
hold off;



return




%------------------------Evaluation of Phase I------------------------------------------------------------------

function dydt = phaseI(t, y, theta)
dydt    = zeros(7,1);

% Parameters
B       = theta(1);              % 1 / person day
c       = theta(2);              % 
p       = theta(3);              % 
N       = theta(4);              % 
r1      = theta(5);              % 1 / days
r2      = theta(6);              % 1 / days
r3      = theta(7);              % 1 / days
r4      = theta(8);              % 1 / days
n       = theta(9);              % 
u       = theta(10);             % 1 / day
h       = theta(11);             % 
alpha   = theta(12);             % 1 / days
v0      = theta(13);             % 
v1      = theta(14);             % 
del     = theta(15);             % 
f       = theta(16);             % 

% Variables
S0      = y(1);
I0_1    = y(2);
I0_2    = y(3);
I0_3    = y(4);
I0_4    = y(5);
Z       = y(6);
D       = y(7);    

% Equations
dydt(1)     = -B * I0_3 * S0;
dydt(2)     =  B * I0_3 * S0 - r1 * I0_1;
dydt(3)     = r1 * I0_1 - r2 * I0_2;
dydt(4)     = r2 * I0_2 - r3 * I0_3; 
dydt(5)     = r3 * I0_3 - r4 * I0_4; 
dydt(6)     = (1 - del) * r4 * (I0_4);
dydt(7)     = del * r4 * (I0_4);

return;


%------------------------Evaluation of Phase II ------------------------------------------------------------------

function dydt = phaseII(t, y, theta)
dydt    = zeros(17,1);


% Parameters
B       = theta(1);              % 1 / person day
c       = theta(2);              % 
p       = theta(3);              % 
N       = theta(4);              % 
r1      = theta(5);              % 1 / days
r2      = theta(6);              % 1 / days
r3      = theta(7);              % 1 / days
r4      = theta(8);              % 1 / days
n       = theta(9);              % 
u       = theta(10);             % 1 / day
h       = theta(11);             % 
alpha   = theta(12);             % 1 / days
v0      = theta(13);             % 
v1      = theta(14);             % 
del     = theta(15);             % 
f       = theta(16);             % 

% Variables
S0      = y(1);
I0_1    = y(2);
I0_2    = y(3);
I0_3    = y(4);
I0_4    = y(5);
Q_0     = y(6);
Q_1     = y(7);
Q_2     = y(8);
Q_3     = y(9);
H       = y(10);
S1      = y(11);
I1_1    = y(12);
I1_2    = y(13);
I1_3    = y(14);
I1_4    = y(15);
Z       = y(16);
D       = y(17);

I3      = I0_3 + Q_3 + I1_3;
Q       = Q_0 + Q_1 + Q_2 + Q_3;
R0      = B * (S0 + Q_0 + S1) / r3;             
K       = (c - p * R0) * r3 * I3 / N;
q2      = (r1 / (r1 + r3 + K)) * (r3 + K) / (r2 + r3 + K);
q3      = ((r1 / (r1 + r3 + K))) * (r2 / (r2 + r3 + K)) * (r3 + K) / (r3 + r3 + K);
q1      = (r3 + K) / (r1 + r3 + K);
L1      = q1 * B * S0 / (r3 + K);
L2      = q2 * B * S0 / (r3 + K);
L3      = q3 * B * S0 / (r3 + K);

if (Q < 1e-6)
    prob = 1;
else
    prob = min([1, n/Q]);
end 


% Equations
dydt(1)     = -B * I3 * S0 - [c - p * R0] * S0 * r3 * I3 / N;
dydt(2)     =  B * I3 * S0 - ((c - p * R0) * (I0_1 / N) + p * L1) * r3 * I3 - r1 * I0_1;
dydt(3)     = r1 * I0_1 -((c - p * R0) * (I0_2 / N) + p * L2) * r3 * I3 - r2 * I0_2;
dydt(4)     = r2 * I0_2 -((c - p * R0) * (I0_3 / N) + p * L3) * r3 * I3 - r3 * I0_3; 
dydt(5)     = r3 * I0_3 - r4 * I0_4; 
dydt(6)     = (c - p * R0) * S0 * r3 * I3 / N - B * I3 * Q_0 - u * Q_0 * prob;
dydt(7)     = B * I3 * Q_0 + ((c - p * R0) * (I0_1 / N) + p * L1) * r3 * I3 - u * Q_1 * prob - r1 * Q_1;
dydt(8)     = r1 * Q_1 + ((c - p * R0) * (I0_2 / N) + p * L2) * r3 * I3 - u * Q_2 * prob - r2 * Q_2;
dydt(9)     = r2 * Q_2 + ((c - p * R0) * (I0_3 / N) + p * L3) * r3 * I3 - u * Q_3 * prob - r3 * Q_3;
dydt(10)    = (1 - f) * h * u * Q_3 * prob - r3 * H - alpha * H;
dydt(11)    = (1 - f) * (1 - v0) * u * Q_0 * prob - B * S1 * I3;
dydt(12)    = B * S1 * I3 + (1 - f) * (1 - v1) * u * Q_1 * prob - r1 * I1_1; 
dydt(13)    = r1 * I1_1 + (1 - f) * u * Q_2 * prob - r2 * I1_2;
dydt(14)    = r2 * I1_2 + (1 - f) * (1 - h) * u * Q_3 * prob + alpha * H - r3 * I1_3;
dydt(15)    = r3 * (I1_3 + Q_3 + H) - r4 * I1_4;
dydt(16)    = (1 - f) * (v0 * Q_0 + v1 * Q_1) * u * prob + (1 - del) * r4 * (I0_4 + I1_4);
dydt(17)    = f * u * Q * prob + del * r4 * (I0_4 + I1_4);

return;