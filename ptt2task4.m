% ASEN 3802 Heat Conduction Lab
% Task 4 - Model II (Updated Initial Temperature Distribution)

clear
clc
close all

%% Geometry
L = 0.1016;
x = 0.0127:0.0127:L;

%% Data files
files = {
    "Aluminum_25V_240mA", "Aluminum_30V_290mA", "Brass_25V_237mA", "Brass_30V_285mA", "Steel_22V_203mA"};

names = { ...
    "Aluminum 25V 240mA", ...
    "Aluminum 30V 290mA", ...
    "Brass 25V 237mA", ...
    "Brass 30V 285mA", ...
    "Steel 22V 203mA"};

%% Material properties (correct from table)
k   = [130 130 115 115 16.2];    % W/(m·K)
rho = [2810 2810 8500 8500 8000];
cp  = [960 960 380 380 500];     % J/(kg·K)

%% Steady-state slopes 
H_exp = [54.907 78.573 104.605 150.712 286.461];

%% Initial slopes 
M_exp = [10 14 18 10 46];

%% Initial intercepts
T0 = [17 16.5 19 19 15];

figure

for f = 1:length(files)

    data = readmatrix(files{f});

    alpha = k(f)/(rho(f)*cp(f));

    H = H_exp(f);
    M = M_exp(f);
    T_initial = T0(f);

    t = linspace(0,max(data(:,1)),500);

    %% Number of Fourier terms
    N = 50;

    T_analytical = zeros(length(x),length(t));

    for i = 1:length(x)

        for j = 1:length(t)

            sum_term = 0;

            for n = 1:N

                lambda_n = (2*n-1)*pi/(2*L);

                % Model II Fourier coefficient
                b_n = (-1)^n*(8*L*(M-H))/((2*n-1)^2*pi^2);
               
                sum_term = sum_term + b_n*sin(lambda_n*x(i))*exp(-lambda_n^2*alpha*t(j));

            end

            % Analytical temperature
            T_analytical(i,j) = T_initial + H*x(i) - sum_term;

        end

    end

    %% Plot results
    subplot(2,3,f)
    hold on

    colors = lines(8);

    for i = 1:8

        % Analytical
        plot(t,T_analytical(i,:), ...
            'Color',colors(i,:), ...
            'LineWidth',2, ...
            'LineStyle','-');

        % Experimental
        plot(data(:,1),data(:,i+1), ...
            'Color',colors(i,:), ...
            'LineWidth',2, ...
            'LineStyle','--');

    end

    xlabel('Time (s)')
    ylabel('Temperature (°C)')
    title(['Model II - ' names{f}])
    grid on

end

%% Legend
h1 = plot(nan,nan,'k-','LineWidth',2);
h2 = plot(nan,nan,'k--','LineWidth',2);

legend([h1 h2],{'Analytical','Experimental'}, ...
       'Position',[0.85 0.1 0.1 0.05])