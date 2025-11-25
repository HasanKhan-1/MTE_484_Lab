%% Figure 2: Inner Loop Plant Output (Theta)
figure(2)
% This line now plots the signal feeding the 'out_cntrl_out' block (Theta)
% ASSUMPTION: 'step_ru*wsol' holds the Theta data. 
theta_response = step_ru*wsol; 
plot(T_outer*(1:K), theta_response);
xlabel('Time [s]');
% Using LaTeX for the Theta symbol in the label
ylabel('$\Theta[k]$ (Inner Loop Plant Output)', 'Interpreter', 'latex'); 
title('Signal from out\_cntrl\_out Scope (Theta)');

% Add t=7s vertical line to the plot
hold on;
jhat = round(7/T_outer); 
t_settle = T_outer * jhat;
line([t_settle t_settle], [min(theta_response) max(theta_response)], ...
     'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);

% --- ADDING THE +/- 2% SETTLING BOUNDS ---
% 1. Determine the Steady State Value for Theta (assumed to be 0)
% If Theta settles to a non-zero value, replace 0 with that value (e.g., theta_ss)
theta_ss = 0; 

% 2. Define the Upper and Lower Bounds (±2% of the Steady State Value)
theta_upper = theta_ss + 0.02; % 0 + 0.02
theta_lower = theta_ss - 0.02; % 0 - 0.02

% 3. Plot the Upper Bound line (starting from t=0 for full context)
line([0 plot_time(end)], [theta_upper theta_upper], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1);

% 4. Plot the Lower Bound line (starting from t=0 for full context)
line([0 plot_time(end)], [theta_lower theta_lower], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1);

hold off;