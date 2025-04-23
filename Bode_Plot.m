clc; 
clear;
close all;

%--------Bode Plot-----------
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
% Plot step response of T(s)
draw_step(T_feedback, 'Closed-Loop System T(s)');

%--------Q5-------- locations of the poles
draw_poles(T_feedback);

%--------Q8-------- Ramp Response
[ess, r_t_out, r_y_out] = draw_ramp(T_feedback, 700+200, 700);

%--------Q9-------- Frequency Response
[Gm, Pm, Wgc, Wpc] = draw_Bode_Plot(G_S*H_S);



%-----------Functions------

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

function [wn, zeta,response_info] = draw_step(sys, sys_name)
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
    disp(['Over shoot MP: ', num2str(100*(response_info.peak_response.value-1)), ...
          '% at t = ', num2str(response_info.peak_response.time), ' sec']);
    % Damping characteristics (for complex poles)
    if ~isreal(response_info.poles)
        [wn, zeta] = damp(sys);
        fprintf('Damping ratio (ζ): %.3f\n', zeta(1));
        fprintf('Natural frequency (ωn): %.3f rad/s\n', wn(1));
    end
    
    disp(['Settling time (2%): ', num2str(response_info.settling_time), ' sec']);
    disp(['Rise time (10-90%): ', num2str(response_info.rise_time), ' sec']);
    disp(['Steady-state value: ', num2str(response_info.steady_state)]);
end

function [poles] = draw_poles(sys)
% DRAW_POLES Plots pole-zero map and returns system poles
%   [poles] = draw_poles(sys)
%
%   Input:
%       sys - Transfer function (tf object) or state-space model
%
%   Output:
%       poles - Array of system poles
%
%   Displays:
%       - Pole-zero plot
%       - Pole locations in command window
%       - Stability information

    % Create figure
    figure;
    
    % Plot pole-zero map
    pzmap(sys);
    title(['Pole-Zero Map of: ' inputname(1)]);
    grid on;
    
    % Get poles
    poles = pole(sys);
    
    % Display poles
    disp(['Poles of ' inputname(1) ':']);
    disp(poles);
    
    
    % Damping characteristics (for complex poles)
    if ~isreal(poles)
        [wn, zeta] = damp(sys);
        fprintf('Damping ratio (ζ): %.3f\n', zeta(1));
        fprintf('Natural frequency (ωn): %.3f rad/s\n', wn(1));
    end

end

function [ess, t_out, y_out] = draw_ramp(sys, t_end, zoom_time)
% DRAW_RAMP Plots ramp response in three subplots
%   [ess, t_out, y_out] = draw_ramp(sys, t_end, zoom_time)
%
%   Inputs:
%       sys - Closed-loop transfer function (tf object)
%       t_end - End time for simulation (default: 100 sec)
%       zoom_time - Time to zoom in (default: 700 sec)
%
%   Outputs:
%       ess - Steady-state error
%       t_out - Time vector
%       y_out - System response vector
%
%   Generates figure with three subplots:
%       1. Ideal ramp input
%       2. System response
%       3. Zoomed comparison at specified time

    % Set defaults if not provided
    if nargin < 2
        t_end = 100;
    end
    if nargin < 3
        zoom_time = 700;
    end

    % Create time vector
    t = 0:0.1:t_end;
    
    %getting the ramp
    ramp = tf(1,[1 0]);
    
    % Get response data
    [y_sys, t_sys] = step(sys.*ramp, t);
    [y_ideal, t_ideal] = step(ramp, t);
    
    % Create figure with three subplots
    figure;
    
    % Subplot 1: Ideal ramp input
    subplot(2,1,1);
    plot(t_ideal, y_ideal, 'b');
    hold on;
    plot(t_sys, y_sys, 'r--');
    title('Ramp Response');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    legend('Ideal', 'System', 'Location', 'northwest');
    grid on;
    hold off;
    
    % Subplot 2: Zoomed comparison
    subplot(2,1,2);
    plot(t_ideal, y_ideal, 'b');
    hold on;
    plot(t_sys, y_sys, 'r--');
    xlim([zoom_time-5 zoom_time+5]);
    title(['Zoomed Comparison at t = ', num2str(zoom_time), ' sec']);
    xlabel('Time (sec)');
    ylabel('Amplitude');
    legend('Ideal', 'System', 'Location', 'northwest');
    grid on;
    hold off;
    
    % Calculate steady-state error (use last 10% of simulation)
    final_idx = round(0.9*length(t_sys)):length(t_sys);
    ess = mean(y_ideal(final_idx) - y_sys(final_idx));
    
    % Display results
    disp(['Steady-state error (ess): ', num2str(ess)]);
    
    % Return output data if requested
    if nargout > 1
        t_out = t_sys;
        y_out = y_sys;
    end
end

function [Gm, Pm, Wgc, Wpc] = draw_Bode_Plot(sys)
% BODE_PLOT Analyzes system stability margins and compares margin() 
%   Bode_Plot(sys)
%
%   Input:
%       sys - Transfer function (tf object or state-space model)
%   Outputs:
%       Gm - Gain margin (dB)
%       Pm - Phase margin (degrees)
%       Wgc - Gain crossover frequency (rad/sec)
%       Wpc - Phase crossover frequency (rad/sec)

    % Create margin plot
    figure;
    margin(sys);
    grid on;
    
    % Get stability margins
    [Gm, Pm, Wgc, Wpc] = margin(sys);
    
    % Display results
    disp(['=== Stability Margins for ' inputname(1) ' ===']);
    disp(['Gain Margin: ', num2str(Gm), ' dB at ', num2str(Wgc), ' rad/s']);
    disp(['Phase Margin: ', num2str(Pm), '° at ', num2str(Wpc), ' rad/s']);
    
end