%  Name and date and other comments


function done = hw3_template()
  clear all;  close all; clc;
  done = 0;
  %%---------------------- Generate Figures 2 a-b ------------------------------

  %% Varied to find balanced promoter strengths giving a bistable toggle network
  du0_v   = linspace(0, 25, 100);
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
  saveas(gcf,'Fig2A.png')
  hold off;

  %% Varied to find an imbalanced promoter strength giving a monstable toggle network
  du0_v   = linspace(0, 25, 100);
  du0_u   = 100 ./ (1 + du0_v .^ 3);

  dv0_u   = du0_u;
  dv0_v   = 5 ./ (1 + dv0_u .^ 1);

  figure;
  plot(du0_u, du0_v);
  hold on;
  plot(dv0_u, dv0_v, 'r');
  title('One Steady-states in a balanced network');
  xlabel('u');
  ylabel('v');
  axis([0 120 0 20]);
  legend('du/dt = 0', 'dv/dt = 0');
  saveas(gcf,'Fig2B.png')
  hold off;
  % ---------------Generating Figures 2 c-d: Varying alpha1 and alpha2 to see bifurgation points------------------------------

  pts     = 70;
  alpha1  = logspace(0.1,3,pts);
  alpha2  = logspace(0.1,3,pts);
  beta    = [3, 2, 1.1];
  gamma   = [3, 2, 1.1];
  map     = zeros(pts,pts);
  figure_num = 0;
  for indx = 1:length(beta)
    figure_num = figure_num + 1;
    for m=1:pts
        for n=1:pts
            du0_v   = linspace(0,50,100);
            du0_u   = alpha1(m) ./ (1 + du0_v .^ beta(indx));
            dv0_u   = du0_u;
            dv0_v   = alpha2(n) ./ (1 + dv0_u .^ gamma(indx));
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
    xlabel('log(alpha 2)');
    ylabel('log(alpha 1)');
    title(['Monostable and Bistable Steady-States.   Beta = Gamma = ', num2str(beta(indx))]);
    saveas(gcf, ['Figure 2C-D Monostable and Bistable Steady-States, b=g=', num2str(beta(indx)),'.png'])
  end

  %---------------------------ODE Figure 5a------------------------------------

  %Define Parameters

  a1  = 156.25;
  a2  = 15.6;
  b   = 2.5;
  g   = 1;
  n   = [2.0015, 0.1, 1, 3];
  K   = 2.9618e-5;
  for i=1:length(n)
    theta   = [a1, a2, b, g, n(i), K, 0];
    IPTG    = logspace(-6, -2, 50);
    tspan   = linspace(0, 17, 50);
    init    = [0, 0];
    initb    = [0, 5];

    for k=1:50
        % Solve ODEs for toggle plasmid, Figure 5a
        theta(7) = IPTG(k);
        [t, y]   = ode23s(@eqns, tspan, init, [], theta);
        u(k)     = y(length(t), 1);
        v(k)     = y(length(t), 2);
    end
    % Normalize to max value if needed such as
    u = u ./ max(u);
    v = v ./ max(v);

    for k=1:50
        % Solve ODEs for toggle plasmid, Figure 5a
        theta(7) = IPTG(k);
        [tb, yb]   = ode23s(@eqns, tspan, initb, [], theta);
        ub(k)     = yb(length(tb), 1);
        vb(k)     = yb(length(tb), 2);

    end
    ub = ub ./ max(ub);
    vb = vb ./ max(vb);

    for k=1:50
        % Solve ODEs for Control, Figure 5a
        theta(7) = IPTG(k);
        [tc, yc]   = ode23s(@eqns_Con, tspan, [0], [], theta);
        vc(k)     = yc(length(tc));
    end
    vc = vc ./ max(vc);
    figure('pos',[10 10 1000 600]);
    semilogx(IPTG, v, 'r');
    hold on;
    semilogx(IPTG, vb, 'k');
    semilogx(IPTG, vc, 'b');
    if i == 1
      title(['Figure 5a: Steady-state gene expression after 17-h induction, gamma', num2str(n(i))]);
      xlabel('[IPTG] (M)')
      ylabel('Normalized GFP expression')
      legend('pTAK117 toggle (low) - stable s.s.','pTAK117 toggle (high) - unstable s.s.','pTAK102 control');
      saveas(gcf, ['Fig5a_g=',num2str(n(i)), '.png'])
    else
      title(['Modified Figure 5a: Steady-state gene expression after 17-h induction, gamma=', num2str(n(i))]);
      xlabel('[IPTG] (M)')
      ylabel('Normalized GFP expression')
      legend('pTAK117 toggle (low) - stable s.s.','pTAK117 toggle (high) - unstable s.s.','pTAK102 control');
      saveas(gcf, ['Fig5a_mod_g=',num2str(n(i)), '.png'])
    end
    hold off;
  end

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
  dydt(1)     = a1 / (1 + v^b) - u;   % finish the equation from the paper
  dydt(2)     = a2 / (1 + (u / (1 + (IPTG / K))^n)^g) - v;   % finish the equation from the paper

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
