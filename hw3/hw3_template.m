%  Name and date and other comments


function done = hw3_template()
clear all;  close all;
done = 0;


%---------------------- Generate Figures 2 a-b ------------------------------

% Varied to find balanced promoter strengths giving a bistable toggle network
du0_v   = linspace(0,25,100);
du0_u   = 100 ./ (1 + du0_v .^ 3);

dv0_u   = du0_u;
dv0_v   = 20 ./ (1 + dv0_u .^ 1);
 
figure;
plot(du0_u, du0_v);
hold on;
plot(dv0_u, dv0_v, 'r');
title('Three Steady-states in a balanced network');
xlabel('u');
ylabel('v');
axis([0 120 0 20]);
legend('du/dt = 0', 'dv/dt = 0');
hold off;


% Varied to find an imbalanced promoter strength giving a monstable toggle network







% ---------------Generating Figures 2 c-d: Varying alpha1 and alpha2 to see bifurgation points------------------------------

pts     = 70;
alpha1  = logspace(0.1,3,pts);
alpha2  = logspace(0.1,3,pts);
beta    = 3;
gamma   = 3;
map     = zeros(pts,pts);

for m=1:pts
    for n=1:pts
        
        du0_v   = linspace(0,50,100);
        du0_u   = alpha1(m) ./ (1 + du0_v .^ beta);
        dv0_u   = du0_u;
        dv0_v   = alpha2(n) ./ (1 + dv0_u .^ gamma);    
        count   = 0;
        
        for p=1:99
            if (((du0_v(p) - dv0_v(p)) * (du0_v(p+1) - dv0_v(p+1))) < 0)
                count = count + 1;
            end
        end
        
        if (count > 1)
            count = 3;
        else
            count = 1;
        end
        
        map(n,m) = count;
    end
end

figure;
contourf(log10(alpha2),log10(alpha1),map);
xlabel('log(alpha 2)'); ylabel('log(alpha 1)'); title('Monostable and Bistable Steady-States.   Beta = Gamma = 3');


% ---------------Vary parameters to see other bifurgation points----------------------------------------

beta    = 2;
gamma   = 2;
map     = zeros(pts,pts);


% Add in here new corresponding calculation and new figure%


beta    = 1.1;
gamma   = 1.1;
map     = zeros(pts,pts);


% Add in here new corresponding calculation and new figure%


%---------------------------ODE Figure 5a------------------------------------
 
%Define Parameters
a1  = ;
a2  = ;
b   = ;
g   = ;
n   = ;
K   = ;
theta   = [a1, a2, b, g, n, K, 0]; 
IPTG    = logspace();
tspan   = linspace();
init    = [ ];

for k=1:

    % Solve ODEs for toggle plasmid, Figure 5a
    theta(7) = IPTG(k);
    [t, y]   = ode23s(@eqns, tspan, init, [], theta);
 
    u(k)     = y();  
    v(k)     = y();
    
end

% Normalize to max value if needed such as
m1  = max(u);


figure;
semilogx();

initb    = [ ];
for k=1:

    % Solve ODEs for toggle plasmid, Figure 5a
    theta(7) = IPTG(k);
    [tb, yb]   = ode23s(@eqns, tspan, initb, [], theta);
  
    
end


hold on;
semilogx();


for k=1:

    % Solve ODEs for Control, Figure 5a
    theta(7) = IPTG(k);
    [tc, yc]   = ode23s(@eqns_Con, tspan, [0], [], theta);
 
    
end




semilogx();
title('')
xlabel('')
ylabel('')
legend('');
hold off;


%--------------------------------------- Your Modification / Perturbation ----------------------------------------------







done = 1;

return;


%-----------------------------------Evaluation of Eqns------------------------------------------------------------------

function dydt = eqns(t, y, theta);
dydt    = zeros(2,1);

% Parameters
a1      = theta(1);             
a2      = theta(2);              
b       = theta(3);              
g       = theta(4);              
n       = theta(5);             
K       = theta(6);   
IPTG    = theta(7);

% Variables
u       = y(1);
v       = y(2);

% Equations
dydt(1)     = ;   % finish the equation from the paper
dydt(2)     = ;   % finish the equation from the paper

return;


%------------------------Evaluation of Eqns for Control ------------------------------------------------------------------

function dydt = eqns_Con(t, y, theta);
dydt    = zeros(1,1);

% Parameters
a1      = theta(1);             
a2      = theta(2);              
b       = theta(3);              
g       = theta(4);              
n       = theta(5);             
K       = theta(6);   
IPTG    = theta(7);

% Variables
v       = y(1);

% Equations
dydt(1)     = a2 / (1 + ((a1/((1 + (IPTG/K))^n))^g)) - v;

return;

