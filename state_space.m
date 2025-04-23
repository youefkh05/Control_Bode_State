clc; 
clear;
close all;

%----------------State Space---------
% Given system matrices
A = [0 1; -6 -5];
B = [0; 1];
C = [1 0];
D = [0];
n = 2;               % System order
sys = ss(A,B,C,D);   % State Space model
x0 = [0; 1];         % Initial condition

% Q2: Transfer function conversion
[num, den] = ss2tf(A,B,C,D);
syms s
TF_Manual = C*inv(s*eye(n)-A)*B + D
TF_builtin = tf(num,den)

% Q3: State transition matrix calculation
% Compute Φ(s) = [sI - A]^-1
Phi_s = inv(s*eye(n) - A);

% Compute Φ(t) by inverse Laplace transform
syms t
Phi_t = ilaplace(Phi_s);

% Verify Φ(0) = I
Phi_0 = subs(Phi_t, t, 0);

% Display results
disp('State transition matrix in s-domain (Φ(s)):');
pretty(Phi_s)

disp('State transition matrix in time domain (Φ(t)):');
pretty(Phi_t)

disp('Verification of Φ(0) = I:');
disp(Phi_0);


% Q4: Verify that Φ̇(t) = AΦ(t)
Phi_dot = diff(Phi_t, t);  % Take time derivative of Φ(t)
A_Phi = A*Phi_t;           % Multiply A with Φ(t)

disp('Time derivative of state transition matrix (Φ̇(t)):');
pretty(Phi_dot)

disp('A*Φ(t):');
pretty(A_Phi)

disp('Verification successful: Φ̇(t) = AΦ(t)');

% Q5 Check Controllability and Observability
% Check Controllability
Co = ctrb(A, B);  % Controllability matrix
rank_Co = rank(Co);
disp('Controllability Matrix:');
disp(Co);
disp(['Rank of Controllability Matrix: ', num2str(rank_Co)]);

if rank_Co == n
    disp('System is Controllable (as expected)');
else
    disp('System is Not Controllable (unexpected for this system)');
end

% Check Observability
Ob = obsv(A, C);  % Observability matrix
rank_Ob = rank(Ob);
disp('Observability Matrix:');
disp(Ob);
disp(['Rank of Observability Matrix: ', num2str(rank_Ob)]);

if rank_Ob == n
    disp('System is Observable (as expected)');
else
    disp('System is Not Observable (unexpected for this system)');
end

% Q6: Unforced (Homogeneous) Response
disp('=== Unforced Response Analysis ===');

% Compute state solution x(t) = Φ(t)*x0
x_t = Phi_t * x0;

disp('Unforced state solution x(t):');
pretty(x_t)

% Compute output solution y(t) = C*x(t) + D*u(t)
% Since u(t)=0 for unforced response:
y_t = C*x_t + D*0;

disp('Unforced output response y(t):');
pretty(y_t)

% Plot the results
t_vals = linspace(0, 5, 500);  % Time vector from 0 to 5 seconds

% Convert symbolic expressions to numeric functions
x1_func = matlabFunction(x_t(1));
x2_func = matlabFunction(x_t(2));
y_func = matlabFunction(y_t);

% Evaluate solutions
x1_vals = arrayfun(x1_func, t_vals);
x2_vals = arrayfun(x2_func, t_vals);
y_vals = arrayfun(y_func, t_vals);

% Plot state responses
figure;
subplot(2,1,1);
plot(t_vals, x1_vals, 'b', 'LineWidth', 2);
hold on;
plot(t_vals, x2_vals, 'r--', 'LineWidth', 2);
title('Unforced State Response');
xlabel('Time (s)');
ylabel('State Values');
legend('x_1(t)', 'x_2(t)');
grid on;

% Plot output response
subplot(2,1,2);
plot(t_vals, y_vals, 'm', 'LineWidth', 2);
title('Unforced Output Response y(t)');
xlabel('Time (s)');
ylabel('Output y(t)');
grid on;

% Compare with MATLAB's built-in initial() function
[~,t_num,x_num] = initial(sys,x0,t_vals(end));
y_num = x_num*C';  % Equivalent to C*x since D=0

% Display symbolic solutions
disp(' ');
disp('Analytic Solutions:');
disp('x1(t) = '); pretty(x_t(1))
disp('x2(t) = '); pretty(x_t(2))
disp('y(t) = '); pretty(y_t)

% Q7: Forced Response Analysis (Unit Step Input)
disp('=== Forced Response Analysis ===');

% Using Frequency Domain Approach
U_s = 1/s;  % Laplace transform of unit step
U_t =ilaplace(U_s);

% Compute forced component in frequency domain
X_forced_s = Phi_s * B * U_s;

% Convert to time domain
x_forced_t = ilaplace(X_forced_s);

% Total solution (homogeneous + forced)
x_total_t = x_t + x_forced_t;

% Output solution
y_total_t = C*x_total_t + D*U_t;  % D*u(t) where u(t)=1 for t>0

disp('Forced state solution (from step input):');
pretty(x_forced_t)

disp('Total state solution (unforced + forced):');
pretty(x_total_t)

disp('Total output solution (unforced + forced):');
pretty(y_total_t)

% Direct evaluation using subs()
x1_vals = double(subs(x_total_t(1), t, t_vals));
x2_vals = double(subs(x_total_t(2), t, t_vals));
y_vals = double(subs(y_total_t, t, t_vals));

% Plot results
figure;

% State responses
subplot(2,1,1);
plot(t_vals, x1_vals, 'b', 'LineWidth', 2);
hold on;
plot(t_vals, x2_vals, 'r--', 'LineWidth', 2);
title('Total State Response (Step Input)');
xlabel('Time (s)');
ylabel('State Values');
legend('Analytic x_1(t)', 'Analytic x_2(t)');
grid on;

% Output response
subplot(2,1,2);
plot(t_vals, y_vals, 'm', 'LineWidth', 2);
hold on;
title('Total Output Response y(t) (Step Input)');
xlabel('Time (s)');
ylabel('Output y(t)');
legend('Analytic y(t)');
grid on;

% Display final steady-state values
ss_x1 = limit(x_total_t(1), t, inf);
ss_x2 = limit(x_total_t(2), t, inf);
ss_y = limit(y_total_t, t, inf);

disp(' ');
disp('Steady-State Values:');
disp(['x1(∞) = ' char(ss_x1)]);
disp(['x2(∞) = ' char(ss_x2)]);
disp(['y(∞) = ' char(ss_y)]);


% Q8: State Feedback Design
disp('=== State Feedback Design ===');

% Original system step response
figure;
step(TF_builtin);
title('Original System Step Response');
grid on;

% Design specifications
zeta_desired = 0.7;      % Desired damping ratio
ts_desired = 1;          % Desired settling time (sec)

% Hand analysis to determine desired poles
wn = 4/(zeta_desired*ts_desired);  % Natural frequency from settling time
sigma = zeta_desired*wn;           % Real part of poles
wd = wn*sqrt(1-zeta_desired^2);    % Imaginary part

% Desired characteristic polynomial
desired_poly = (s + sigma + 1i*wd)*(s + sigma - 1i*wd);
desired_poly = expand(desired_poly);

% Convert to numerical polynomial
desired_coeffs = sym2poly(desired_poly);

% Hand calculation of K matrix
% Characteristic polynomial of A-BK: s^2 + (5+K2)s + (6+K1)
% Compare with desired polynomial: s^2 + 2*zeta*wn*s + wn^2

K1 = desired_coeffs(3) - 6;  % From constant term
K2 = desired_coeffs(2) - 5;  % From s term
K = [K1 K2]

disp('Desired closed-loop poles:');
disp([-sigma+1i*wd, -sigma-1i*wd]);

disp('Feedback gain matrix K:');
disp(K);

% Verification
Ac = A - B*K;
[num_2, denum_2] = ss2tf(Ac,B,C,D);
TF_state_feedback = tf(num_2, denum_2);

% Step response analysis
figure;
step_info = stepinfo(TF_state_feedback);
step(TF_state_feedback);
title('System with State Feedback');
grid on;

disp('Closed-loop system performance:');
disp(['Settling Time: ', num2str(step_info.SettlingTime), ' sec']);
disp(['Overshoot: ', num2str(step_info.Overshoot), '%']);
