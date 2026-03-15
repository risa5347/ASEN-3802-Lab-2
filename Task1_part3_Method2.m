% ASEN 3802 Heat Conduction Lab
% Model III - Adjusted Thermal Diffusivity using RMS

clear;clc;close all;

%% Geometry
L = 0.1016;
x = 0.0127:0.0127:L;

%% Files
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

%% Material Properties
k   = [960 960 380 380 500];
rho = [2810 2810 8500 8500 8000];
cp  = [130 130 115 115 16.2];

%% Experimental Slopes
H_exp = [54.907 78.573 104.605 150.712 286.461];

%% Initial Temperatures
T0 = [17 16.5 19 19 15];

%% Diffusivity scaling range
scale = 0.3:0.01:1.5;

alpha_original = k./(rho.*cp);
alpha_adjusted = zeros(size(alpha_original));
min_rms = zeros(length(files),1);

figure

%% Loop through datasets
for f = 1:length(files)

    data = readmatrix(files{f});
    t_data = data(:,1);

    H = H_exp(f);
    T_initial = T0(f);

    rms_values = zeros(length(scale),1);

    %% Search for best alpha
    for s = 1:length(scale)

        alpha = alpha_original(f)*scale(s);

        N = 10;
        T_model = zeros(length(t_data),8);

        for tc = 1:8
            for j = 1:length(t_data)

                sum_term = 0;

                for n = 1:N
                    lambda_n = (2*n-1)*pi/(2*L);
                    b_n = (-1)^n*(8*H*L)/((2*n-1)^2*pi^2);

                    sum_term = sum_term + ...
                        b_n*sin(lambda_n*x(tc))*exp(-lambda_n^2*alpha*t_data(j));
                end

                T_model(j,tc) = T_initial + H*x(tc) + sum_term;

            end
        end

        error = data(:,2:9) - T_model;

        rms_values(s) = sqrt(mean(error(:).^2,'omitnan'));

    end

    %% Best diffusivity
    [min_rms(f),idx] = min(rms_values);
    alpha_adjusted(f) = alpha_original(f)*scale(idx);

    %% Analytical solution with best alpha
    alpha = alpha_adjusted(f);

    t = linspace(0,max(data(:,1)),500);
    N = 10;

    T_analytical = zeros(length(x),length(t));

    for i = 1:length(x)
        for j = 1:length(t)

            sum_term = 0;

            for n = 1:N
                lambda_n = (2*n-1)*pi/(2*L);
                b_n = (-1)^n*(8*H*L)/((2*n-1)^2*pi^2);

                sum_term = sum_term + ...
                    b_n*sin(lambda_n*x(i))*exp(-lambda_n^2*alpha*t(j));
            end

            T_analytical(i,j) = T_initial + H*x(i) + sum_term;

        end
    end

    %% Plot
    subplot(2,3,f)
    hold on
    colors = lines(8);

    for i= 1:8

        plot(t,T_analytical(i,:),'Color',colors(i,:),'LineWidth',2)

        plot(data(:,1),data(:,i+1),'--','Color',colors(i,:),'LineWidth',2)

    end

    xlabel('Time (s)')
    ylabel('Temperature (°C)')
    title(['Model III - ' names{f}])
    grid on

end

sgtitle('Model III: Adjusted Thermal Diffusivity')

%% Global Legend
h1 = plot(nan,nan,'k-','LineWidth',2);
h2 = plot(nan,nan,'k--','LineWidth',2);

legend([h1 h2],{'Analytical','Experimental'},'Position',[0.85 0.1 0.1 0.05])

%% Compute RMS percent error

T_exp_mean = zeros(length(files),1);

for f = 1:length(files)

    data = readmatrix(files{f});
    T_exp_mean(f) = mean(data(:,2:9),'all','omitnan');

end

rms_percent= (min_rms ./ T_exp_mean) * 100;

%% Display table

T = table(names', alpha_original', alpha_adjusted', min_rms, rms_percent, ...
    'VariableNames', {'Material','Alpha_Original','Alpha_Adjusted','Min_RMS_Error','RMS_Error_Percent'});

disp(T)