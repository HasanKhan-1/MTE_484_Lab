%% Plotting the solution
figure(1)
plot_time = T_outer*(1:K);
output_response = step_ry*[xsol;xhatsol];
plot(plot_time, output_response, 'b', 'LineWidth', 1.5);
hold on;
% --- CONSTRAINT LINES ---
% 1. Steady State Reference Line
% The final steady-state value is constrained to be -steadyState(1,:)*[x;xhat]
y_ss = -steadyState(1,:)*[xsol;xhatsol]; 
line([0 plot_time(end)], [y_ss y_ss], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

% 2. Settling Time Lines
% t = 7s corresponds to jhat
jhat = round(7/T_outer); 
t_settle = T_outer * jhat;
% Settling Time Vertical Line (t=7s)
line([t_settle t_settle], [min(output_response) max(output_response)], ...
     'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);
text(t_settle + 0.5, max(output_response) * 0.9, 't = 7s', 'Color', 'r');

% Settling Time Upper and Lower Bounds (±2% of Steady State)
y_upper = 1.02 * y_ss;
y_lower = 0.98 * y_ss;
line([t_settle plot_time(end)], [y_upper y_upper], 'Color', 'm', 'LineStyle', '-.', 'LineWidth', 1);
line([t_settle plot_time(end)], [y_lower y_lower], 'Color', 'm', 'LineStyle', '-.', 'LineWidth', 1);

% 3. Overshoot Constraint Line (1.45 * |Steady State|)
% Overshoot constraint is max(y) <= 1.45 * |y_ss|
y_overshoot_limit = 1.45 * abs(y_ss);
line([0 plot_time(end)], [y_overshoot_limit y_overshoot_limit], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 1);
% --------------------------
xlabel('Time [s]');
ylabel('y[k]');
title('System Output Step Response with Constraints');
grid on;
hold off;

%% Figure 2: Inner Loop Plant Output (Theta)
figure(2)
% This line now plots the signal feeding the 'out_cntrl_out' block (Theta)
% ASSUMPTION: 'step_ru*wsol' holds the Theta data. 
% If your Theta data is in a different variaXble (e.g., 'theta_sol'), replace 'step_ru*wsol'
theta_response = step_ru*wsol; 
plot(T_outer*(1:K), theta_response);
xlabel('Time [s]');
% Using LaTeX for the Theta symbol in the label
ylabel('$\Theta[k]$ (Inner Loop Plant Output)', 'Interpreter', 'latex'); 
title('Signal from out\_cntrl\_out Scope (Theta)');

% Add t=7s vertical line to the plot
hold on;
% Calculate t_settle again just in case, though it's already done above
jhat = round(7/T_outer); 
t_settle = T_outer * jhat;
line([t_settle t_settle], [min(theta_response) max(theta_response)], ...
     'Color', 'r', 'LineStyle', ':', 'LineWidth', 1);
hold off;

% use log scale for heat map?
% ... (rest of the code remains the same)z