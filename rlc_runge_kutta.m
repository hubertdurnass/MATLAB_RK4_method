clc; clear;

%% @file rlc_runge_kutta.m
%  @brief Solves the second-order differential equation for an RLC circuit using RK4.
%  @author Hubert Durnas

%% Circuit Parameters
% @brief Defines the component values for the RLC circuit.
L = 0.1; % @var L Inductance in Henrys (H)
C = 0.01; % @var C Capacitance in Farads (F)
R = 5 * 2 * sqrt(L / C); % @var R Resistance in Ohms (Ω) - underdamped case

%% Square Wave Parameters
% @brief Parameters for the input square wave signal.
A = rand(1) * 10; % @var A Amplitude of the square wave
freq = 0.5; % @var freq Frequency of the square wave in Hz

%% Initial Conditions
% @brief Defines the starting conditions for voltage and current.
u0 = 0; % @var u0 Initial voltage [V]
v0 = 0; % @var v0 Initial current [A]
Y0 = [u0; v0]; % @var Y0 Initial state vector

%% Time Step and Simulation Time
h = 0.001; % @var h Time step [s]
T = 5.0; % @var T Total simulation time [s]
N = round(T / h); % @var N Number of steps

%% Initialize Results
% @brief Arrays for storing simulation data.
t = linspace(0, T, N + 1); % @var t Time vector
solution = zeros(2, N + 1); % @var solution Solution matrix for voltage and current
solution(:, 1) = Y0;

U = zeros(1, N + 1); % @var U Input signal array

%% Fourth-Order Runge-Kutta Method
% @brief Implements the RK4 numerical integration method.
for i = 1:N
    % @brief Generates input signal as a square wave.
    U(i) = A * square(2 * pi * freq * t(i));

    k1 = f(t(i), solution(:, i), R, L, C, U(i));
    k2 = f(t(i) + 0.5 * h, solution(:, i) + 0.5 * h * k1, R, L, C, U(i));
    k3 = f(t(i) + 0.5 * h, solution(:, i) + 0.5 * h * k2, R, L, C, U(i));
    k4 = f(t(i) + h, solution(:, i) + h * k3, R, L, C, U(i));

    solution(:, i + 1) = solution(:, i) + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end

%% Test using ODE45 Solver
% @brief Compares RK4 results with MATLAB's built-in ode45 solver.
ode45Function = @(t, Y) f(t, Y, R, L, C, A * square(2 * pi * freq * t));
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10); % @var options Higher precision settings
[t_ode45, solution_ode45] = ode45(ode45Function, [0, T], Y0, options);

%% Output Results
% @brief Extracts voltage output values for visualization.
Uwy = solution(1, :); % @var Uwy Voltage output from RK4
Uwy_test = solution_ode45(:, 1); % @var Uwy_test Voltage output from ode45 solver

%% Plot Results
% @brief Generates a plot comparing the results of RK4 and ode45.
figure;
plot(t, U, 'g-', 'LineWidth', 2, 'DisplayName', 'Input signal U');
hold on;
plot(t, Uwy, 'b-', 'LineWidth', 2, 'DisplayName', 'Runge-Kutta method');
plot(t_ode45, Uwy_test, 'r--', 'LineWidth', 2, 'DisplayName', 'ODE45 solver');
xlabel('Time [s]');
ylabel('Amplitude');
title('Output Voltage Uwy(t)');
legend('show'); grid on;
