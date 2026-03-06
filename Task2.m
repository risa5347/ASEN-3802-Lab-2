% ASEN 3802 Heat Conduction Lab
% Analytical vs Experimental Temperature vs Time
clear
clc 
close all

%% Material Properties
Brass = readmatrix("Brass_25V_237mA");
k = 115;
rho = 8500;
cp = 380;
alpha = k/(rho*cp);

%% Geometry
L = 0.1016;

% Thermocouple locations (m)
x = 0.0127 : 0.0127 : L;  % Th1 -> Th8

%% Steady State Parameters from Part 1
H_an = 101.68;
T0 = 19;

%% Time Vector
t = linspace(0,max(Brass(:,1)),500);

%% Number of Modes
N = 10;
modes_array = 0:N;

%% Analytical Temperature Matrix

%% Analytical Temperature Matrix
T_analytical = zeros(length(x), length(t));

for i = 1:length(x)          % loop over thermocouples / spatial positions
    for j = 1:length(t)      % loop over time
        sum_term = 0;
        for n = 1:N
            lambda_n = (2*n - 1)*pi / (2*L);
            b_n = (8 * H_an * L) / ((2*n - 1)^2 * pi^2) * (-1)^n;
            sum_term = sum_term + b_n * sin(lambda_n * x(i)) * exp(-lambda_n^2 * alpha * t(j));
        end
        T_analytical(i,j) = T0 + H_an * x(i) + sum_term;   % steady + transient
    end
end
%% Plot
figure
hold on

% Define colors for 8 thermocouples
colors = lines(8);

for i = 1:8
    % Analytical (solid)
    plot(t, T_analytical(i,:), 'Color', colors(i,:), 'LineWidth', 2, 'LineStyle', '-');
    
    % Experimental (dashed)
    plot(Brass(:,1), Brass(:,i+1), 'Color', colors(i,:), 'LineWidth', 2, 'LineStyle', '--');
end

xlabel('Time (s)')
ylabel('Temperature (°C)')
title('Model IA: Analytical vs Experimental Temperature vs Time')
grid on

% Custom legend for line styles
h1 = plot(nan, nan, 'k-', 'LineWidth', 2);   % Solid black (analytical)
h2 = plot(nan, nan, 'k--', 'LineWidth', 2);  % Dashed black (experimental)
legend([h1, h2], {'Analytical', 'Experimental'}, 'Location', 'northwest')