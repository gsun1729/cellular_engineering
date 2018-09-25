 % Gordon Sun
 % 20180924
 % HW2 MAPK Kinase from "Ultrasensitivity in mitogen activated protein kinase cascade"

 % Sections of this code come from Upinder S. Bhalla and NCBS

function n = SunGordon_HW2()
	% clear all;  close all; clc;
	% Using default run time from 0 to 100
	tspan = [0, 30];
	% Initial conditions of 14 unbuffered molecules + 10 enzyme-substrate complexes
	y0 = [0.0003; 0.1; 1.2; 0; 0; 1.2; 0; 0; 0.003; 0; 0.0003; 0.0003; 0.12; 0.12; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];

	% Molecule names
	% 1	E2
	% 2	INPUT(E1)
	% 3	MAPK
	% 4	MAPK-P
	% 5	MAPK-PP
	% 6	MAPKK
	% 7	MAPKK-P
	% 8	MAPKK-PP
	% 9	MAPKKK
	% 10	MAPKKK*
	% 11	MAPKKPase
	% 12	MAPKKPase[1]
	% 13	MAPKPase
	% 14	MAPKPase[1]
	% 15	E2_complex
	% 16	INPUT(E1)_complex
	% 17	MAPKK-PP_complex
	% 18	MAPKK-PP_MAPKK-PP[1]_complex
	% 19	MAPKKK*_complex
	% 20	MAPKKK*_MAPKKK*[1]_complex
	% 21	MAPKKPase_complex
	% 22	MAPKKPase[1]_MAPKKPase_complex
	% 23	MAPKPase_complex
	% 24	MAPKPase[1]_MAPKPase_complex

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   Define # of s.s. points and stimulus range here
	%   A different stimulus range can be used for each enzyme (see paper)

	% # of s.s. evaluations
	pts = 100;
	% A different stimulus range will be used for each enzyme
	stim = zeros(pts, 3);
	% MAPKhas a narrow range due to ultrasensitivity
	MAPK = logspace(-6, -3, pts)';
	% MAPKK has a medium range
	MAPKK = logspace(-6, -2, pts)';
	% MAPKKK has a broad stimulus range
	MAPKKK = logspace(-6, 0, pts)';

	input_vals = [MAPK, MAPKK, MAPKKK];
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   Get your s.s. enzyme concentrations from the ODE solver here
	%   Likely want a loop / a few loops for solving over stimulus range
	%   Stimulus is changed by changing the corresponding initial condition

	% values of key parameters at s.s.
	v = zeros(pts, 3);

	for index = 1:pts
		% MAPK-P Index 5
		y0(2)   = input_vals(index, 1);
		[t, y]  = ode23s(@f, tspan, y0);
		col     = length(y(:, 1));
		v(index, 1) = y(col, 5);
		% MAPKK-PP INdex 8
		y0(2)   = input_vals(index, 2);
		[t, y]  = ode23s(@f, tspan, y0);
		col     = length(y(:, 1));
		v(index, 2) = y(col, 8);
		% MAPKKK* index 10
		y0(2)   = input_vals(index, 3);
		[t, y]  = ode23s(@f, tspan, y0);
		col     = length(y(:, 1));
		v(index, 3) = y(col, 10);
    end

 	%   Recommendation:  Normalize enzyme concentrations by max concentrations
	normalized_dat = zeros(pts, 3);
    [pts, enzymes] =size(normalized_dat);

	for n_index = 1:enzymes
		normalized_dat(:, n_index) = v(:, n_index) ./ v(pts, n_index);
    end
	% Recommendation:  Use interpolation to find EC10 and EC90 points
	% = spline(Y, X, [a b]);
	EC_pts(1, :) = spline(normalized_dat(:, 1), MAPK, [0.10, 0.50, 0.90]);
	EC_pts(2, :) = spline(normalized_dat(:, 2), MAPKK, [0.10, 0.50, 0.90]);
	EC_pts(3, :) = spline(normalized_dat(:, 3), MAPKKK, [0.10, 0.50, 0.90]);

	% Hill Coefficient Calculation
	n_h = log(81) ./ log(EC_pts(:, 3) ./ EC_pts(:, 1));

	% Plots (Make sure all figures are labeled with legends and axes titles)
	% Concentrations vs. time
	Font = 12;
    figure;
	plot(t, [ y(:, 5), y(:, 8), y(:, 10), y(:, 3), y(:, 6), y(:, 9)]);
    xlabel('Time (min)', 'FontSize', Font);
    ylabel('Conc. (uM)', 'FontSize', Font);
    title('Species Conc. vs Time for [E1]_0 = 0.1uM', 'FontSize', Font);
    legend('MAPK-PP','MAPKK-PP','MAPKKK*','MAPK','MAPKK','MAPKKK');
	saveas(gcf,'speciesConc_time.png')
	% Plots for Figure2B
	figure;

	semilogx(MAPK, normalized_dat(:, 1));
	hold on
	semilogx(MAPKK, normalized_dat(:, 2),'r');
	semilogx(MAPKKK, normalized_dat(:, 3),'g');
	plot([10^-6 1],[0.1 0.1], 'black');
	plot([10^-6 1],[0.9 0.9], 'black');
	legend('MAPK', 'MAPKK', 'MAPKKK');
	title('Figure 2B', 'FontSize', Font);
	xlabel('Input Stimulus E1_{tot} (nM)', 'FontSize', Font);
	ylabel('Predicted S.S. Kinase Activity', 'FontSize', Font);
	saveas(gcf,'Fig2B.png')
	% Plot for Hill Coefficient Calculation (Figure 2A)
	figure;
	hold on;
	plot(MAPK ./ EC_pts(1, 2), normalized_dat(:, 1), 'b');
	plot(MAPKK ./ EC_pts(2, 2), normalized_dat(:, 2), 'r');
	plot(MAPKKK ./ EC_pts(3, 2), normalized_dat(:, 3), 'g');
	title('Figure 2A', 'FontSize', Font);
	xlabel('Input Stimulus E1_{tot} (nM)', 'FontSize', Font);
	ylabel('Predicted S.S. Kinase Activity', 'FontSize', Font);
	axis([0 6 0 1]);

	% Evaluate Hill Eqn Curves of Comparable Steepness to Plot
	% These plots are just overlaid onto the Figure 2A plots (can be black
	% curves)

	%Recommendation:
	for indx = 1:enzymes
		hill_curve_fit(:, indx) = Hill_equation([n_h(indx), EC_pts(indx, 2)], input_vals(:, indx));
	end
	plot(MAPK ./ EC_pts(1, 2), hill_curve_fit(:, 1), 'k-.');
	plot(MAPKK ./ EC_pts(2, 2), hill_curve_fit(:, 2), 'k-.');
	plot(MAPKKK ./ EC_pts(3, 2), hill_curve_fit(:, 3), 'k-.');
	legend('MAPK','MAPKK','MAPKKK','Hill Eqn Curves');
	saveas(gcf,'Fig2A.png')
return;

%------------------------Evaluation function------------------------------------------------------------------
function dydt = f(t, y)
	dydt = [
	- 1000 * y(10) * y(1)+ 300 * y(15)
	- 1000 * y(9) * y(2)+ 300 * y(16)
	- 1000 * y(3)* y(8)+ 150 * y(17)+ 150 * y(23)
	- 1000 * y(4)* y(8)+ 150 * y(18)- 1000 * y(4)* y(13)+ 150 * y(23)+ 150 * y(17)+ 150 * y(24)
	- 1000 * y(5)* y(14)+ 150 * y(24)+ 150 * y(18)
	- 1000 * y(6)* y(10)+ 150 * y(19)+ 150 * y(21)
	- 1000 * y(7)* y(10)+ 150 * y(20)- 1000 * y(7)* y(11)+ 150 * y(21)+ 150 * y(19)+ 150 * y(22)
	- 1000 * y(8)* y(12)+ 150 * y(22)+ 150 * y(20)- 1000 * y(3) * y(8)+ 300 * y(17)- 1000 * y(4) * y(8)+ 300 * y(18)
	- 1000 * y(9)* y(2)+ 150 * y(16)+ 150 * y(15)
	- 1000 * y(10)* y(1)+ 150 * y(15)+ 150 * y(16)- 1000 * y(6) * y(10)+ 300 * y(19)- 1000 * y(7) * y(10)+ 300 * y(20)
	- 1000 * y(7) * y(11)+ 300 * y(21)
	- 1000 * y(8) * y(12)+ 300 * y(22)
	- 1000 * y(4) * y(13)+ 300 * y(23)
	- 1000 * y(5) * y(14)+ 300 * y(24)
	1000 * y(10) * y(1)- 300 * y(15)
	1000 * y(9) * y(2)- 300 * y(16)
	1000 * y(3) * y(8)- 300 * y(17)
	1000 * y(4) * y(8)- 300 * y(18)
	1000 * y(6) * y(10)- 300 * y(19)
	1000 * y(7) * y(10)- 300 * y(20)
	1000 * y(7) * y(11)- 300 * y(21)
	1000 * y(8) * y(12)- 300 * y(22)
	1000 * y(4) * y(13)- 300 * y(23)
	1000 * y(5) * y(14)- 300 * y(24)
	];

return;


%------------------Function that describes the curve that data points are fit to (Optional) --------------
function yhat = Hill_equation(beta, E_free)
	n = beta(1)
	K_d = beta(2);
	yhat = ((E_free).^n) ./ (K_d.^n + (E_free).^n);
return;
