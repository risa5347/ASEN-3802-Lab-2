clear;
clc;
close all;

Brass25 = readmatrix("Steel_22V_203mA");
Brass25_Time = Brass25(:,1);

Brass = Brass25(:,2:9); % all 8 datasets
K_Brass = 150;
A_Brass = ((pi/4)*1)*0.00064516; %m^2
I_Brass = 237*10^(-3); %A
V_Brass = 25; %V
Q_dot = I_Brass*V_Brass;
H_Brass = Q_dot/(K_Brass*A_Brass); %K/m

% figure;
% hold on;
% 
% for i = 1:8
%     y = Brass(:,i);
% 
%     plot(Brass25_Time, y, '-');
% 
%     p = polyfit(Brass25_Time, y, 1);
%     y_fit = polyval(p, Brass25_Time);
% 
%     plot(Brass25_Time, y_fit, '-');
% end
% 
% xlabel('Time');
% ylabel('Value');
% title('Brass Data with Line of Best Fit');
% grid on;
% 
% figure;
% hold on;
% 
% for i = 1:8
%     y = Brass(:,i);
% 
%     plot(Brass25_Time, y, '-');
% 
%     p = polyfit(Brass25_Time, y, 1);
%     y_fit = polyval(p, Brass25_Time);
% 
%     plot(Brass25_Time, y_fit, '-');
% end
% 
% xlabel('Time');
% ylabel('Value');
% title('Brass Data with Line of Best Fit');
% grid on;
% 
% hold off;

figure;
hold on;

startTime = Brass25_Time > 10000;

err_val = 0.5; % <-- CHANGE THIS (your measurement uncertainty)

for i = 1:8
    x = Brass25_Time(startTime);
    y = Brass(startTime,i);

    err = err_val * ones(size(y));  % constant error

    % Plot with error bars
    errorbar(x, y, err, '-');

    % Best fit
    p = polyfit(x, y, 1);
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);

    plot(x_fit, y_fit, '-');
end

% xlabel('Time');
% ylabel('Value');
% title('Brass, 25V, 237mA Data with Error Bars & Fit');
% grid on;
% 
% figure;
% hold on;
% 
% % --- 1. Get steady-state temperatures (average after 1500) ---
% startTime = Brass25_Time > 1500;
% T_ss = mean(Brass(startTime,:), 1);   % 1x8 vector
% 
% % --- 2. Define sensor positions (YOU MUST EDIT THESE) ---
% % Example: equally spaced (replace with real values from your lab)
% x_pos = linspace(0, .5, 8);  % meters (CHANGE THIS)
% 
% % --- 3. Error bars (same idea as before) ---
% err_val = 0.2;
% err = err_val * ones(size(T_ss));
% 
% % --- 4. Plot experimental data ---
% errorbar(x_pos, T_ss, err, 'o', 'LineWidth', 1.5);
% 
% % --- 5. Experimental best-fit line ---
% p_exp = polyfit(x_pos, T_ss, 1);
% x_fit = linspace(min(x_pos), max(x_pos), 100);
% y_exp_fit = polyval(p_exp, x_fit);
% 
% plot(x_fit, y_exp_fit, '-', 'LineWidth', 2);
% 
% % --- 6. Analytical line ---
% % T = H*x + C → need an intercept (use first point as reference)
% T0 = T_ss(1);
% y_analytical = H_Brass * x_fit + T0;
% 
% plot(x_fit, y_analytical, '--', 'LineWidth', 2);
% 
% % --- Labels ---
% xlabel('Position (m)');
% ylabel('Temperature');
% title('Steady-State Temperature Distribution (Experimental vs Analytical for Steel 22V 203mA)');
% legend('Experimental (with error)', 'Experimental Fit', 'Analytical', 'Location', 'best');
% grid on;
% 

figure;
hold on;

% --- Initial temperatures (first recorded time step) ---
T_initial = Brass(1,:);   % 1x8 vector

% --- Same sensor positions ---
x_pos = linspace(0, .5, 8);  % use your real positions if different

% --- Error bars (optional, same idea as before) ---
err_val = 0.05;
err = err_val * ones(size(T_initial));

% --- Plot ---
errorbar(x_pos, T_initial, err, 'o', 'LineWidth', 1.5);

% Optional: connect points for visualization
plot(x_pos, T_initial, '-', 'LineWidth', 1.5);

xlabel('Position (m)');
ylabel('Temperature');
title('Initial Temperature Distribution (All Steel 22V 203mA Thermocouples)');
grid on;
hold off;