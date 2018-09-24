 % Gordon Sun
 % 20180924
 % HW2 MAPK Kinase from "Ultrasensitivity in mitogen activated protein kinase cascade"

 % Sections of this code come from Upinder S. Bhalla and NCBS

function n = HW2_Template()
	clear all;  close all; clc;
	% Using default run time from 0 to 100
	tspan=[0, 100];
	% Initial conditions of 14 unbuffered molecules + 10 enzyme-substrate complexes
	y0 = [0.0003; 0.1; 1.2; 0; 0; 1.2; 0; 0; 0.003; 0; 0.0003; 0.0003; 0.12; 0.12; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

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
	pts = 50;         
	% A different stimulus range will be used for each enzyme
	stim = zeros(pts, 3);
	% MAPKhas a narrow range due to ultrasensitivity
	stim(:,1) = logspace(-7, -4, pts)';
	% MAPKK has a medium range
	stim(:,2) = logspace(-7, -4, pts)';
	% MAPKKK has a broad stimulus range
	stim(:,3) = logspace(-7, -1, pts)';
	% values of key parameters at s.s.
	v = zeros(pts,3);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   Get your s.s. enzyme concentrations from the ODE solver here
	%   Likely want a loop / a few loops for solving over stimulus range
	%   Stimulus is changed by changing the corresponding initial condition

		[t, y]  = ode23s(@f, tspan, y0);    % Main solver to call

	  
		% One example for MAPK:
		% disp(' ');  disp('Calculating MAPK'); 
	for p = 1:pts
		y0(2)   = stim(p, 1);
		[t, y]  = ode23s(@f, tspan, y0);
		col     = length(y(:,1));
		v(p,1)  = y(col,5);              % 5	    MAPK-PP  .

		y0(2)   = stim(p, 1);
		[t, y]  = ode23s(@f, tspan, y0);
		col     = length(y(:,1));
		v(p,1)  = y(col,5);              % 5	    MAPK-PP  

		y0(2)   = stim(p, 1);
		[t, y]  = ode23s(@f, tspan, y0);
		col     = length(y(:,1));
		v(p,1)  = y(col,5);              % 5	    MAPK-PP  
	 end 
		
		% Similarly think of MAPKK and MAPKKK
		
		
		
		
		
		
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
		


	%   Recommendation:  Normalize enzyme concentrations by max concentrations




	% Recommendation:  Use interpolation to find EC10 and EC90 points
	% = spline(Y, X, [a b]);



	% Use interpolation to find EC50 points


	% Hill Coefficient Calculation 
	% n =                    % The equation is found in the lecture notes



	% Plots (Make sure all figures are labeled with legends and axes titles)
	% Concentrations vs. time
	figure;
	plot(t, y);


	% Plots for Figure2B
	% figure;
	% semilogx();
	% hold on;
	% semilogx(,'r');
	% semilogx(,'g');
	% plot([10^-6 0.1],[0.1 0.1], 'black');
	% plot([10^-6 0.1],[0.9 0.9], 'black');
	% legend();
	% title();
	% xlabel();
	% ylabel();

	% Plot for Hill Coefficient Calculation (Figure 2A)
	% hold off;
	% figure;
	% plot();
	% hold on;
	% plot(,'r');
	% plot(,'g');
	% title();
	% xlabel();
	% ylabel();
	% axis([]);

	% Evaluate Hill Eqn Curves of Comparable Steepness to Plot
	% These plots are just overlaid onto the Figure 2A plots (can be black
	% curves)

	%Recommendation:
	%for 1:3
	%
	%end
	% plot();
	% legend();
end;

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