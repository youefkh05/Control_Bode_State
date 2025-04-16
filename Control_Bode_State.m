clc; 
clear;
close all;

%--------Q1--------  Define G(s) and H(s)
% Define the open-loop transfer function G(s)
num_G = 1;
den_G = [1 1 0];    % s(s+1) = s^2 + s

% Define the feedback transfer function H(s)
num_H = [1];
den_H = [1];        % Unity Feedback
[G_S, H_S] = create_system(num_G, den_G, num_H, den_H)

%--------Q2-------- Step Response of G(s) (Open-Loop)
% Plot step response of G(s)
draw_step(G_S, 'Open-Loop System G(s)');

%--------Q3-------- Closed-Loop Analysis
disp('Closed-Loop TF using feedback():');
T_feedback = feedback(G_S, H_S)

disp('Closed-Loop TF using manual formula (G/(1+GH)):');
T_manual =   (1 / (1 + G_S * H_S)) * G_S;  % Equivalent to T(s) = G/(1+GH)
T_manual = minreal(T_manual)               % Cancel common terms

%--------Q4-------- Step Response of T(s) (Closed-Loop)
% Plot step response of G(s)
draw_step(T_feedback, 'Closed-Loop System T(s)');


%-----Functions------

function [G_S, H_S] = create_system(num_G, den_G, num_H, den_H)
% CREATE_SYSTEM Creates open-loop and feedback transfer functions
%   [G_S, H_S] = create_system(num_G, den_G, num_H, den_H)
%   
%   Inputs:
%       num_G - Numerator coefficients of G(s)
%       den_G - Denominator coefficients of G(s)
%       num_H - Numerator coefficients of H(s) (default: 1)
%       den_H - Denominator coefficients of H(s) (default: 1)
%
%   Outputs:
%       G_S - Open-loop transfer function
%       H_S - Feedback transfer function

    % Set default unity feedback if not specified
    if nargin < 3
        num_H = 1;
        den_H = 1;
    end
    
    % Create transfer functions
    G_S = tf(num_G, den_G)
    H_S = tf(num_H, den_H)
end

function [response_info] = draw_step(sys, sys_name)
% DRAW_STEP Plots step response and returns key performance metrics
%   [response_info] = draw_step(sys, sys_name)
%
%   Inputs:
%       sys - Transfer function (tf object)
%       sys_name - Name of the system for title (string)
%
%   Outputs:
%       response_info - Structure containing:
%           .poles - System poles
%           .stability - Stability classification
%           .peak_response - Peak response value and time
%           .settling_time - Time to settle within 2% of final value
%           .rise_time - 10-90% rise time
%           .steady_state - Final steady-state value
%       Figure with step response

    % Create figure
    figure;
    
    % Get step response data
    [y, t] = step(sys);
    
    % Plot step response
    step(sys);
    title(['Step Response of ', sys_name]);
    grid on;
    
    % Calculate response characteristics
    response_info = struct();
    response_info.poles = pole(sys);
    
    % Stability determination
    if all(real(response_info.poles) < 0)
        response_info.stability = 'stable (all poles in LHP)';
    elseif any(real(response_info.poles) > 0)
        response_info.stability = 'unstable (at least one pole in RHP)';
    else
        response_info.stability = 'marginally stable (poles on imaginary axis)';
    end
    
    % Peak response (overshoot)
    [response_info.peak_response.value, peak_idx] = max(y);
    response_info.peak_response.time = t(peak_idx);
    
    % Steady-state value (last 10% of response)
    steady_state_val = mean(y(end-round(length(y)*0.1):end));
    response_info.steady_state = steady_state_val;
    
    % Settling time (within 2% of steady-state)
    settled_idx = find(abs(y - steady_state_val) > 0.02*steady_state_val, 1, 'last');
    if isempty(settled_idx)
        response_info.settling_time = 0;
    else
        response_info.settling_time = t(settled_idx);
    end
    
    % Rise time (10% to 90% of steady-state)
    rise_start = find(y >= 0.1*steady_state_val, 1);
    rise_end = find(y >= 0.9*steady_state_val, 1);
    if ~isempty(rise_start) && ~isempty(rise_end)
        response_info.rise_time = t(rise_end) - t(rise_start);
    else
        response_info.rise_time = NaN;
    end
    
    % Display results in command window
    disp(['System: ', sys_name]);
    disp(['Poles: ', num2str(response_info.poles')]);
    disp(['Stability: ', response_info.stability]);
    disp(['Peak response: ', num2str(response_info.peak_response.value), ...
          ' at t = ', num2str(response_info.peak_response.time), ' sec']);
    disp(['Settling time (2%): ', num2str(response_info.settling_time), ' sec']);
    disp(['Rise time (10-90%): ', num2str(response_info.rise_time), ' sec']);
    disp(['Steady-state value: ', num2str(response_info.steady_state)]);
end