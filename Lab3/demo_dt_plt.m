%{
1.Start program in Arduino
2.Go to file Connections -> File Capture -> Start
3. Save data to where this script is saved
4. Once the data is completed save it with good_name.txt
5. Replace txt files with name in script
6. Run it !
%}


% --- Configuration ---
data_file = 'data3.txt';
Ts = 0.027646; % Sampling interval (s) from your Arduino code
% --- Data Loading ---
fprintf('Loading data from %s...\n', data_file);
% Attempt to load the data file. The data must be in CSV format (e.g., 0.1000,0.0912).
try
   % Read the CSV data. Data should be loaded into a matrix [y_ref, y_current].
   data = readmatrix(data_file);
catch
   error('Could not read the data file. Please ensure "%s" is accessible and correctly formatted (CSV).', data_file);
end
% Separate columns
y_ref = data(:, 1);
y_current = data(:, 2);
N = length(y_current); % Number of samples
fprintf('Found %d data points.\n', N);
% --- Time Vector Calculation (The "Adder Block") ---
% The time vector is calculated as: time = [0, 1*Ts, 2*Ts, ..., (N-1)*Ts]
time_vector = (0:(N-1)) * Ts;
% --- Plotting ---
figure('Name', 'Control System Step Response Analysis');
hold on;
grid on;
grid minor;
% Plot the Reference (Input)
plot(time_vector, y_ref, 'b--', 'LineWidth', 2, 'DisplayName', 'Reference (y\_ref)');
% Plot the Current Position (Output)
plot(time_vector, y_current, 'r-', 'LineWidth', 2, 'DisplayName', 'Measured Output (y\_current)');
% --- Customization and Labels ---
title(sprintf('Sampled-Data System Step Response (Ts = %g s)', Ts));
xlabel('Time (seconds)');
ylabel('Ball Position (m)');
legend('show', 'Location', 'southeast');
set(gca, 'FontSize', 12); % Adjust font size for better readability
hold off;
fprintf('Plotting complete. Check the generated figure window.\n');
