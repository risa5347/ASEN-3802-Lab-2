% ASEN 3802 Heat Conduction Lab
% Analytical vs Experimental Temperature vs Time
clear
clc
close all

% Geometry
L = 0.1016;

% Thermocouple locations (m)
x = 0.0127 : 0.0127 : L;  % Th1 -> Th8

% Files
files = { ...
    "Aluminum_25V_240mA", ...
    "Aluminum_30V_290mA", ...
    "Brass_25V_237mA", ...
    "Brass_30V_285mA", ...
    "Steel_22V_203mA" };

names = { ...
    "Aluminum 25V 240mA", ...
    "Aluminum 30V 290mA", ...
    "Brass 25V 237mA", ...
    "Brass 30V 285mA", ...
    "Steel 22V 203mA" };

% Material Properties
k   = [960 960 380 380 500];      
rho = [2810 2810 8500 8500 8000];
cp  = [130 130 115 115 16.2];

% Steady State Parameters
H_an = [54.907 78.573 104.605 150.712 286.461];
%^ Is now used as H_exp after part 2 task 3
T0   = [17 16.5 19 19 15];

% Create Figure Once
figure

% Loop Through Files
for f = 1:length(files)

    data = readmatrix(files{f});

    alpha = k(f)/(rho(f)*cp(f));
    H = H_an(f);
    T_initial = T0(f);

    % Time Vector
    t = linspace(0,max(data(:,1)),500);

    % Number of Modes
    N = 10;

    % Analytical Temperature Matrix
    T_analytical = zeros(length(x), length(t));

    for i = 1:length(x)
        for j = 1:length(t)

            sum_term = 0;

            for n = 1:N
                lambda_n = (2*n - 1)*pi/(2*L);
                b_n = (8*H*L)/((2*n-1)^2*pi^2) * (-1)^n;

                sum_term = sum_term + ...
                    b_n*sin(lambda_n*x(i))*exp(-lambda_n^2*alpha*t(j));
            end

            T_analytical(i,j) = T_initial + H*x(i) + sum_term;

        end
    end

    % Subplot
    subplot(2,3,f)
    hold on

    colors = lines(8);

    for i = 1:8

        plot(t, T_analytical(i,:), ...
            'Color', colors(i,:), 'LineWidth',2,'LineStyle','-');

        plot(data(:,1), data(:,i+1), ...
            'Color', colors(i,:), 'LineWidth',2,'LineStyle','--');

    end

    xlabel('Time (s)')
    ylabel('Temperature (°C)')
    title(names{f})
    sgtitle('IA')
    grid on

end

% Global Legend
h1 = plot(nan,nan,'k-','LineWidth',2);
h2 = plot(nan,nan,'k--','LineWidth',2);

legend([h1 h2],{'Analytical','Experimental'},'Position',[0.85 0.1 0.1 0.05])