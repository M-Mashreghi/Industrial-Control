clc
clear all

%%
k1=0.8
b2=0.8
k2=0.28
b1=0.28
m1=1
m2=1
syms s

X2_to_X1 =(b1*s+k1)/(m2*s^2  + (b2+b1)*s +(k2+k1))

X1_to_U = (m2*s^2  + (b2+b1)*s +(k2+k1))/((m2*s^2  + (b2+b1)*s +(k2+k1))*(m1*s^2+b1*s+k1)-(b1*s+k1)^2)

X2_to_U = X2_to_X1 * X1_to_U
% Simplify X2_to_U
X2_to_U_simplified = simplify(X2_to_U);
display(X2_to_U_simplified)
%% Q3

% Obtain numerator and denominator of X2_to_U_simplified
[numerator, denominator] = numden(X2_to_U_simplified);

% Display numerator and denominator
disp('Numerator of X2_to_U_simplified:');
disp(numerator);
disp('Denominator of X2_to_U_simplified:');
disp(denominator);

% Convert symbolic coefficients to numeric values
coefficients_numerator = fliplr(double(coeffs(numerator, s)));
coefficients_denominator = fliplr(double(coeffs(denominator, s)));

% Display numeric coefficients
disp('Numeric coefficients of numerator:');
disp(fliplr(coefficients_numerator));
disp('Numeric coefficients of denominator:');
disp(fliplr(coefficients_denominator));

% Create transfer function
g = tf((coefficients_numerator), (coefficients_denominator));
display(g);

% Display step response information
stepinfo(g);
display(stepinfo(g))
%%
step(g)
sys_info = stepinfo(g);
rise_time = sys_info.RiseTime;   
T_s = rise_time/10;  
display(T_s)
%%
z = tf('z',T_s);
G_z = c2d(g,T_s,'tustin')
[numerator, denominator] = tfdata(G_z, 'v');
figure()
hold on
step(G_z)
step(g)
title('Discretized vs Continious System')
legend('Discrete','Continious')
hold off
%% T

% Define the number of sinusoidal components
num_components = 6;
G_i = randi([1, 7], 1, 6);% Random amplitudes between 0 and 7
display(G_i)
omega_i = randi([1, 7], 1, 6); % Random angular frequencies between 0 and 7
display(omega_i)
%%

% Create a time vector for simulation
t = 0:T_s:50; % Adjust the time vector as needed

% Initialize the input signal
u = ones(size(t));

% Generate the sum of sinusoidal inputs
for i = 1:num_components
    u = u + G_i(i) * sin(omega_i(i) * t);
end
disp(u)


figure;
plot(t, u);
xlabel('Time');
ylabel('Composite Input Signal');
title('Sum of Sinusoidal Inputs');
grid on;



figure()
[y,t_out] = lsim(G_z,u,t);
% Plot the system response
figure;
plot(t_out, y);
xlabel('Time');
ylabel('System Response');
title('System Response to Sinusoidal Input');
grid on;
%% Se

[parameter,regres_matrix] = least_square(y,u,t);
display(parameter)
y_estimated = regres_matrix*parameter;
0.5*sum((y-y_estimated).^2)


% % Extract numerator and denominator coefficients
% numerator_coeffs = parameter(5:9)';
% denominator_coeffs = [1,parameter(1:4)'];
% 
% % Create transfer function G_seZ in discrete time
% G_seZ = tf(numerator_coeffs,denominator_coeffs, T_s);
% G_se = d2c(G_seZ)
% display(G_se)
%% je

for j = 1:1
    
    % Create a time vector for simulation
    t = 0:T_s:50; % Adjust the time vector as needed
    % Initialize the input signal
    u = ones(size(t))*2;
    % Generate the sum of sinusoidal input
    figure;
    plot(t, u);
    xlabel('Time');
    ylabel('Composite Input Signal');
    title('Sum of Sinusoidal Inputs');
    grid on;

    figure()
    [y,t_out] = lsim(G_z,u,t);
    % Plot the system response
    figure;
    plot(t_out, y);
    xlabel('Time');
    ylabel('System Response');
    title('System Response to Sinusoidal Input');
    grid on;
    [parameter,regres_matrix] = least_square(y,u,t);
    display(parameter)
    y_estimated = regres_matrix*parameter;
    disp("loss")
    disp(0.5*sum((y-y_estimated).^2))
    
end


for j = 1:5
    
    % Create a time vector for simulation
    t = 0:T_s:50; % Adjust the time vector as needed

    % Initialize the input signal
    u = ones(size(t));

    % Generate the sum of sinusoidal inputs
    for i = 1:j
         u = u + G_i(i) * sin(omega_i(i) * t);
    end
    
    figure;
    plot(t, u);
    xlabel('Time');
    title(['Sum of Sinusoidal Inputs for j = ', num2str(j)]);
    title('Sum of Sinusoidal Inputs');
    grid on;

    figure()
    [y,t_out] = lsim(G_z,u,t);
    % Plot the system response
    figure;
    plot(t_out, y);
    xlabel('Time');
    ylabel('System Response');
    title(['System Response to Sinusoidal Input for j = ', num2str(j)]);
    grid on;

    [parameter,regres_matrix] = least_square(y,u,t);
    disp(j+1)
    display(parameter)
    y_estimated = regres_matrix*parameter;
    disp("loss")
    disp(0.5*sum((y-y_estimated).^2))

end

%% #che

% Create a time vector for simulation
t = 0:T_s:50; % Adjust the time vector as needed

% Initialize the input signal
u = ones(size(t));

% Generate the sum of sinusoidal inputs
for i = 1:num_components
    u = u + G_i(i) * sin(omega_i(i) * t);
end



figure;
plot(t, u);
xlabel('Time');
ylabel('Composite Input Signal');
title('Sum of Sinusoidal Inputs');
grid on;

figure()
[y,t_out] = lsim(G_z,u,t);
% Plot the system response
figure;
plot(t_out, y);
xlabel('Time');
ylabel('System Response');
title('System Response to Sinusoidal Input with noise');
grid on;

y_noisy = y+randn(131,1);
[parameter,regres_matrix] = least_square(y_noisy,u,t);
regres_matrix*parameter
0.5*sum((y_noisy-y_estimated).^2)
disp(parameter)
%% he

% Create a time vector for simulation
t = 0:T_s:50; % Adjust the time vector as needed

% Initialize the input signal
u = ones(size(t));

% Generate the sum of sinusoidal inputs
for i = 1:num_components
    u = u + G_i(i) * sin(omega_i(i) * t);
end
[y,~] = lsim(G_z,u,t);
[~,parameters] = r_least_square(y,u,t);
[numerator , denominator] = tfdata(G_z, 'v');
actual_params = [denominator(2:end),numerator];

figure()
for i=1:4
    hold on
    subplot(2,2,i)
    plot(parameters(i,:))
    yline(actual_params(i),'-.','Color','r','LineWidth',1)
    legend("Actual Value","Estimated Value")
    title(['numerator',num2str(i-1)])
    hold off
end
figure()
for i=1:5
    hold on
    subplot(2,3,i)
    plot(parameters(i+4,:))
    yline(actual_params(i+4),'-.','Color','r','LineWidth',1)
    legend("Actual Value","Estimated Value")
    title(['denominator',num2str(i-1)])
    hold off
end

%% khe

% Create a time vector for simulation
t = 0:T_s:50; % Adjust the time vector as needed

% Initialize the input signal
u = ones(size(t));

% Generate the sum of sinusoidal inputs
for i = 1:num_components
    u = u + G_i(i) * sin(omega_i(i) * t);
end
[y,~] = lsim(G_z,u,t);
y_noisy = y+randn(131,1);
[~,parameters] = r_least_square(y_noisy,u,t);
[numerator , denominator] = tfdata(G_z, 'v');
actual_params = [denominator(2:end),numerator];

figure()
for i=1:4
    hold on
    subplot(2,2,i)
    plot(parameters(i,:))
    yline(actual_params(i),'-.','Color','r','LineWidth',1)
    legend("Actual Value","Estimated Value")
    title(['numerator',num2str(i-1)])
    hold off
end
figure()
for i=1:5
    hold on
    subplot(2,3,i)
    plot(parameters(i+4,:))
    yline(actual_params(i+4),'-.','Color','r','LineWidth',1)
    legend("Actual Value","Estimated Value")
    title(['denominator',num2str(i-1)])
    hold off
end