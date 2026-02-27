% ASEN 3802 Heat Conduction Lab - Part 2 Task 1
% Analytical Transient Solution Convergence Plot

clear; clc; close all;

%% 1. Material and Geometric Parameters
% Aluminum 7075-T651 Properties
k = 130;       % Thermal Conductivity [W/(m*K)]
rho = 2810;    % Density [kg/m^3]
cp = 960;      % Specific Heat Capacity [J/(kg*K)]
alpha = k / (rho * cp); % Thermal Diffusivity [m^2/s]

% Geometry Calculations (Converting inches to meters)
% According to the manual: x0 is 1.375 in to the left of Th1.
% Th1 to Th8 consists of 7 spaces of 0.5 in = 3.5 in.
% The heater begins 1.0 in past Th8.
x_Th8 = (1.375 + 3.5) * 0.0254; % Position of Th8 [m]
L = x_Th8 + (1.0 * 0.0254);     % Total length of rod from x0 [m]

%% 2. Part 1 Calculated Parameters
% IMPORTANT: Verify these match your exact values from Part 1 Task 1
H_an = 94.88; % Analytical steady-state slope for Al 25V [deg C/m]
T0 = 19.0;    % Steady-state cold end temp [deg C] (Approx from your Fig 2.2)

%% 3. Time and Modes Setup
times = [1, 1000]; % Times to evaluate [s]
max_modes = 10;    % Evaluate from 0 to 10 modes
modes_array = 0:max_modes;

% Initialize temperature arrays for the two evaluation times
T_t1 = zeros(1, length(modes_array));
T_t2 = zeros(1, length(modes_array));

%% 4. Calculate Temperatures
% Steady-state baseline (which represents n = 0)
T_steady = T0 + H_an * x_Th8;

for i = 1:length(modes_array)
    num_modes = modes_array(i);
    
    sum_t1 = 0;
    sum_t2 = 0;
    
    % Sum the transient terms up to current num_modes
    for n = 1:num_modes
        % Eigenvalue lambda_n
        lambda_n = ((2*n - 1)*pi) / (2*L);
        
        % Derived b_n coefficient from your Part 2 Task 1 derivation
        b_n = (8 * H_an * L) / (((2*n - 1)^2) * pi^2) * (-1)^n;
        
        % Spatial sine term
        spatial_term = sin(lambda_n * x_Th8);
        
        % Transient exponential terms
        transient_term_t1 = exp(-lambda_n^2 * alpha * times(1));
        transient_term_t2 = exp(-lambda_n^2 * alpha * times(2));
        
        % Accumulate summations
        sum_t1 = sum_t1 + (b_n * spatial_term * transient_term_t1);
        sum_t2 = sum_t2 + (b_n * spatial_term * transient_term_t2);
    end
    
    % Total temperature = Steady State + Transient Sum
    T_t1(i) = T_steady + sum_t1;
    T_t2(i) = T_steady + sum_t2;
end

%% 5. Plotting
figure('Name', 'Convergence of Analytical Transient Solution', 'Color', 'w');

% Plot data for t = 1s and t = 1000s
plot(modes_array, T_t1, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 't = 1 s');
hold on;
plot(modes_array, T_t2, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 't = 1000 s');

% Formatting to meet grading requirements
title('Convergence of Analytical Solution at Th_8 (Aluminum 25V)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Number of Modes (n)', 'FontSize', 11);
ylabel('Temperature (^{\circ}C)', 'FontSize', 11);
grid on;
legend('Location', 'best', 'FontSize', 10);
xticks(0:max_modes);

%% 6. Fourier Number Calculation (For Discussion Section)
Fo_t1 = (alpha * times(1)) / (L^2);
Fo_t2 = (alpha * times(2)) / (L^2);

fprintf('--- Fourier Numbers ---\n');
fprintf('Fourier Number at t = 1s: %.4e\n', Fo_t1);
fprintf('Fourier Number at t = 1000s: %.4f\n', Fo_t2);