% Autonomous Vehicles Laboratory. P.I.: PhD Jesús Alberto Sandoval Galarza.
% at Instituto Tecnológico de La Paz (ITLP).

% Undergraduate student: Leonardo Corral Trigueros. May, 2026

% Experiment: Replicating the experiment of the IDA-PBC controller with unbounded torques for the Inertia Wheel Pendulum from: 
% "Control of the Inertia Wheel Pendulum by Bounded Torques" [1], derived from
% "Stabilization of a class of underactuated mechanical systems via interconnection and damping assignment" [2].  

% [1] doi: https://doi.org/10.1016/j.ifacol.2020.12.1749
% [2] doi: 10.1109/TAC.2002.800770


function IDA_PBC_IWP_nbt()
    % Parameters from Section V
    m3 = 10;
    m11 = 0.1;
    m22 = 0.2;
    u_max = 45; % Reference upper limit for comparison
    
    % Derived/Chosen parameters.
    a1 = 1; 
    a2 = -1.5; 
    a3 = 6;
    
    gamma1 = 30; 
    gamma2 = 4.5; 
    %k2 = 0.0266;
    
    % Controller Gains
    kp = 3.75;
    kv = 10;
    
    % Additional derived constants for Equation (16)/(15)
    delta = a1*a3 - a2^2;
    k3 = -(a1 + a2) / delta;
    k4 = -(a2 + a3) / (a1 + a2);

    % Pack parameters into a structure for the ODE function
    params = struct('m3', m3, 'm11', m11, 'm22', m22, ...
                    'gamma1', gamma1, 'gamma2', gamma2, ...
                    'k3', k3, 'k4', k4, 'kp', kp, 'kv', kv);

    % Initial configurations from Section V
    % Rest initial state: [q1(0); q2(0); p1(0); p2(0)] = [3.14; 0; 0; 0]
    x0 = [3.14; 0; 0; 0]; 
    tspan = [0, 40]; % Time horizon: 40 seconds

    % Run simulation
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, x] = ode45(@(t, x) system_dynamics(t, x, params), tspan, x0, options);

    % Compute the control torque post-simulation for plotting
    u = zeros(length(t), 1);
    for i = 1:length(t)
        q1 = x(i, 1);
        q2 = x(i, 2);
        p1 = x(i, 3);
        p2 = x(i, 4);
        
        % Original Control Law (Equation 16) / unbounded torque
        u(i) = gamma1*sin(q1) + kp*(q2 + gamma2*q1) + kv*k3*(p2 + k4*p1);
    end

    % --- Plotting Results ---
    
    % Figure 3 Replication: Joint Positions q1 and q2
    figure('Position', [100, 100, 600, 400]);
    plot(t, x(:, 1), 'k-', 'LineWidth', 1.5); hold on;
    plot(t, x(:, 2), 'k--', 'LineWidth', 1.2);
    grid on;
    xlabel('t [sec]');
    ylabel('[rad]');
    legend('q_1', 'q_2', 'Location', 'best');
    title('IDA-PBC. Figure 3: Joint Positions q_1 and q_2 (Non-Bounded Control Law)');
    ylim([-2, 3.5]);

    % Figure 5 Replication: Applied Non-Bounded Torque u
    figure('Position', [750, 100, 600, 400]);
    plot(t, u, 'k-', 'LineWidth', 1.5); hold on;
    % Plot the u_max threshold boundary for visualization
    plot(tspan, [u_max, u_max], 'r:', 'LineWidth', 1.2);
    plot(tspan, [-u_max, -u_max], 'r:', 'LineWidth', 1.2);
    grid on;
    xlabel('t [sec]');
    ylabel('[Nm]');
    title('IDA-PBC. Figure 5: Applied Non-Bounded Torque u');
    ylim([-55, 55]);
    
    % Display peak torque value achieved
    fprintf('Peak absolute torque achieved: %.2f Nm\n', max(abs(u)));
end


% Desired closed loop acquisition:
function dxdt = system_dynamics(~, x, params)
    % Extract states
    q1 = x(1);
    q2 = x(2);
    p1 = x(3);
    p2 = x(4);

    % Extract parameters
    m3  = params.m3;
    m11 = params.m11;
    m22 = params.m22;
    gamma1 = params.gamma1;
    gamma2 = params.gamma2;
    kp = params.kp;
    kv = params.kv;
    k3 = params.k3;
    k4 = params.k4;

    % Original Control Law: Equation (16)
    u = gamma1*sin(q1) + kp*(q2 + gamma2*q1) + kv*k3*(p2 + k4*p1);

    % Hamiltonian Equations of Motion: Equation (5)
    dq1_dt = p1 / m11;
    dq2_dt = p2 / m22;
    dp1_dt = m3 * sin(q1) - u;
    dp2_dt = u;

    dxdt = [dq1_dt; dq2_dt; dp1_dt; dp2_dt];
end