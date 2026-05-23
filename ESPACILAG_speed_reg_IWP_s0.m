% Autonomous Vehicles Laboratory. P.I.: PhD Jesús Alberto Sandoval Galarza.
% at Instituto Tecnológico de La Paz (ITLP).


% Undergraduate student: Leonardo Corral Trigueros. May, 2026

% Experiment: Comparing the next controller for the Inertia Wheel Pendulum with the same
% conditions (\omega = 0 rad/s) as the controller from: 
    % "Control of the Inertia Wheel Pendulum by Bounded Torques" [1], derived from
    % "Stabilization of a class of underactuated mechanical systems via interconnection and damping assignment" [2].  
    
    % [1] doi: https://doi.org/10.1016/j.ifacol.2020.12.1749
    % [2] doi: 10.1109/TAC.2002.800770

function ESPACILAG_speed_reg_IWP_s0()
    %Controller from "A Speed Regulator for A Torque-Driven Inertia Wheel Pendulum" [3].
    %[3] doi: https://doi.org/10.1016/j.ifacol.2020.12.1749

    % Physical Parameters (Section 2.1 & 3)
    I1 = 0.1;   % Moment of inertia of pendulum 
    I2 = 0.2;   % Moment of inertia of disk 
    m3 = 10;    % Gravity/mass parameter 
    
    a1 = I1 + I2; 
    a2 = I2;       
    a3 = I2;     
    M = [a1, a2; a2, a3]; 
    detM = a1*a3 - a2^2;  
    
    % Controller Gains (Section 3)
    kp = 0.5;   
    kv = 0.15;  
    d1 = 5;     
    d2 = 2;     
    d3 = 1;     
    % Desired constant speed r
    %r_des = 5;  % ESPACILAG speed regulation paper [3]
    r_des = 0; % Comparison to SPACILAG controller based on IDA-PBC equivalent [1]

    Ma = [d1, d2; d2, d3]; 
    detMa = d1*d3 - d2^2;  
    gamma2 = (d1 - d2)/(d3 - d2); 
    
    % Pack parameters into a structure for the ODE function
    params = struct('M', M, 'm3', m3, 'Ma', Ma, 'detM', detM, 'detMa', detMa, ...
                    'kp', kp, 'kv', kv, 'gamma2', gamma2, 'r_des', r_des, ...
                    'a1', a1, 'a2', a2, 'a3', a3, 'd1', d1, 'd2', d2, 'd3', d3); 
    
    % --- Simulation Setup ---
    % Initial conditions: [q1, q2, p1, p2] -> Pendulum at 85.9 deg, rest [3] 
    %x0 = [deg2rad(85.9); 0; 0; 0]; %original x0 ESPACILAG speed regulation paper
    x0 = [3.14; 0; 0; 0]; % test x0 IDA-PBC condition [1]

    tspan = [0 12]; % Time interval 
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6); 
    
    % Solve using ODE45 (Passing the custom dynamics function)
    [t, x] = ode45(@(t, x) dynamics(t, x, params), tspan, x0, options);
    
    % --- Post-Processing Control Inputs ---
    u_computed = zeros(length(t), 1);
    for i = 1:length(t)
        u_computed(i) = compute_espacilag_tau(t(i), x(i, :)', params);
    end
    
    % =========================================================================
    % FIGURE 2 REPLICATION (Modified: q1 in Radians, Wheel Speed in rad/s)
    % =========================================================================
    figure('Color', 'w', 'Name', 'Figure 2 Replication');
    
    % Subplot 1: Pendulum Position q1 
    subplot(2, 1, 1);
    plot(t, x(:,1), 'b', 'LineWidth', 1.2); 
    hold on;
    plot(t, zeros(size(t)), 'g--', 'LineWidth', 1); 
    ylabel('q_1 [rad]');  
    title('q_1 [rad]');   
    grid on;
    %ylim([-100 100]); % original q1 in deg, r_des= 5 [3]
    ylim([-2 4]); %test q1 in rad, r_des = 0 [1]    
    xlim([0 12]);        
    
    % Subplot 2: Wheel Speed dq2/dt 
    subplot(2, 1, 2);
    dq = (M \ x(:, 3:4)')'; 
    q2_dot = dq(:, 2);
    
    plot(t, q2_dot, 'b', 'LineWidth', 1.2); 
    hold on;
    plot(t, r_des*ones(size(t)), 'r', 'LineWidth', 1); 
    plot(t, zeros(size(t)), 'g--', 'LineWidth', 1);    
    ylabel('$$\dot{q}_2$$ [rad/s]', 'Interpreter', 'latex');
    title('$$\dot{q}_2$$ [rad/s]', 'Interpreter', 'latex');   
    xlabel('Time [s]'); 
    grid on;
    %ylim([-10 25)]; % original \dot{q2} in rad/s, r_des = 5 [3]
    ylim([-30 35]); % test \dot{q2} in rad/s, r_des = 0 [1]      
    xlim([0 12]);
    
    % =========================================================================
    % FIGURE 5 REPLICATION (Control Torque Profile)
    % =========================================================================
    %Plotting from [1]
    figure('Color', 'w', 'Position', [680, 100, 550, 400], 'Name', 'Figure 5 Replication');
    plot(t, u_computed, 'b', 'LineWidth', 1.2); hold on;
    plot(t, zeros(size(t)), 'g--', 'LineWidth', 1);
    grid on; 
    xlabel('Time [s]'); 
    ylabel('Control Torque u [Nm]');
    title('ESPACILAG (q, p): Replication of Figure 5.');
    ylim([-15, 15]);
    xlim([0 12]);
end

% Modular Torque Computation Function
function u = compute_espacilag_tau(t, x, params)
    q = x(1:2);
    p = x(3:4);
    
    a1 = params.a1;
    a2 = params.a2;
    a3 = params.a3;
    d1 = params.d1;
    d2 = params.d2;
    d3 = params.d3;
    M = params.M;
    Ma = params.Ma;
    m3 = params.m3;
    detMa = params.detMa;
    gamma2 = params.gamma2;
    kp = params.kp;
    kv = params.kv;
    detM = params.detM;
    r = params.r_des;
    
    % Error Coordinates qa 
    % K is diag(a1, a2), Section 2.5 
    qa = [a1 * q(1); a2 * (q(2) - r*t)]; 
    
    % Calculate pa 
    q_dot = M \ p;
    qa_dot = [a1 * q_dot(1); a2 * (q_dot(2) - r)];
    pa = Ma * qa_dot;
    
    % Control Law components 
    % Potential Shaping u_es 
    dUa_dqa1 = -(m3 * detMa / (d3 - d2)) * sin(qa(1)/a1) - gamma2 * kp * (qa(2) - gamma2 * qa(1));
    dUa_dqa2 = kp * (qa(2) - gamma2 * qa(1));
    
    ues = (1/detMa) * ((-(a2/a1)*d3 + (a3/a2)*d2)*dUa_dqa1 + ((a2/a1)*d2 - (a3/a2)*d1)*dUa_dqa2);
    
    % Damping Injection u_di 
    udi = -(kv * a1 * a2 / detM) * (-pa(1) + pa(2));
    
    % Final control torque u acting on the wheel
    u = ues + udi;
end

% Desired closed loop acquisition:
function dxdt = dynamics(t, x, params)
    q = x(1:2);
    p = x(3:4);
    
    M = params.M;
    m3 = params.m3;
    
    % Get control torque u from the isolated controller function
    u = compute_espacilag_tau(t, x, params);
    
    % State derivatives calculations
    dq_dt = M \ p;                  % Evaluates generalized velocities
    dp_dt = [m3 * sin(q(1)); u];    % Evaluates generalized momentum
    
    dxdt = [dq_dt; dp_dt];
end