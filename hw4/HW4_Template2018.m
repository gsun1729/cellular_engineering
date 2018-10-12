

function y0 = HW4_Template2018()
clear all;  close all;

%---------------------------------------------------SET UP--------------------------------------------------------------

%Set Parameters
B       =              % 1 / person day
c       =                % 
p       =               % 
N       =              % 
r1      =               % 1 / days
r2      =               % 1 / days
r3      =               % 1 / days
r4      =              % 1 / days
n       =             % 
u       =                % 1 / day (TV)
h       =               % 
alpha   =             % 1 / days
v0      =             % 
v1      =            % 
del     =              % 
f       =              % 
I_t0    = 
 
%Pack Parameters
theta = [B, c, p, N, r1, r2, r3, r4, n, u, h, alpha, v0, v1, del, f]; 


% Set initial conditions for Phase I


% Solve ODEs for Phase I


[t, y]  = ode23s( @    ,    ,     ,   , theta );



%------------------------------------------------Trace Vaccination----------------------------------------------------------------

% Set initial conditions for Phase II TV


% Solve ODEs for Phase II TV


[t2, y2]  = ode23s( @  ,    ,     ,  , theta);

% Keep track of infected, queued, and dead

L1   = length(t);
L2   = length(t2);

YourArray = zeros(L1+L2,1);

for m=1:L1
   YourArray(m) = t(m);
   
end
for n=1:L2
    YourArray(n+L1) = t2(n)+5; 
    
end



%------------------------------------------------Mass Vaccination-------------------------------------------------------------------

% Set initial conditions for Phase II MV


theta(10) = 200;

%Solve ODEs for Phase II MV

[t3, y3]  = ode23s( @  ,    ,     ,  , theta);

% Keep track of infected, queued, and dead

L3       = length(t3);

YourArray2 = zeros(L1+L3,1);

for m=1:L1
   YourArray2(m) = t(m);
   
end
for n=1:L3
    YourArray2(n+L1) = t3(n)+5; 
    
end



%-------------------------------------------------------------------------Plots------------------------------------------------

figure;
[AX,H1,H2] = plotyy( YourTime1, YourY1, YourTime2, YourY2);
set(get(AX(1),'Ylabel'),'String','Number in Queue');
set(get(AX(2),'Ylabel'),'String','Number Infected');
set(H1,'LineStyle',':');
set(H2,'LineStyle','-');
title('MyTitle');
xlabel('YourTime');
legend(H1, 'YourY1');
legend(H2, 'YourY2');




%------------------------------------------------Death Rate Comparison------------------------------------------------------------    
    
% CDC policy (CDC guidelines are modeled by switching from TV to MV 28 days after the start of TV (e.g. on day 33 in the base case))


% Set initial conditions 

% Solve

[t4, y4]  = ode23s( @  ,    ,     ,  , theta);


% Set initial conditions 

% Solve

[t5, y5]  = ode23s( @  ,    ,     ,  , theta);


% Combine results 



% Plot TV vs MV vs CDC Policy

figure;
plot( , );
hold on;
plot(   ,    , 'r');
plot(   ,    , 'g');
title('Your Title for Small Pox Deaths');
xlabel('  ');
ylabel('   ');
legend('TV', 'MV', 'CDC Policy');
hold off;




%------------------------Sensitivity Analysis----------------------------

% Fraction of population initially infected

tspan   = [0 5];
tspan2  = [0 400];
dead    = zeros(50,1);
dead2   = zeros(50,1);
frac    = zeros(50,1);
for f=1:50
    frac(f) = % such as f divided by something 
    i0      = % such as frac(f) multiplied by something
    init    = [N i0 0 0 0 0 0];
    [t, y]  = ode23s(@phaseI, tspan, init, [], theta);
    
    init2   = [ ];
    theta(10) = 50;
    [t2, y2] = ode23s(@phaseII, tspan2, init2, [], theta);
    
    dead(f) = y2( , 17);    
    
    init4   = [ ];
    theta(10) = 200;
    [t4, y4]  = ode23s(@phaseII, tspan2, init4, [], theta);
   
    dead2(f)  = y4( , 17);    
    
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








return;


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
