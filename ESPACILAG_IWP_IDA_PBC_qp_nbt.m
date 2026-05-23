% Autonomous Vehicles Laboratory. P.I.: PhD Jesús Alberto Sandoval Galarza.
% at Instituto Tecnológico de La Paz (ITLP).


% Undergraduate student: Leonardo Corral Trigueros. May, 2026

% In order to know the acquisition of the desired closed-loop proposed by ESPACILAG methodology. We recommend checking:
% A note on the energy-shaping control applied to mechanical systems: "A pendulum example [1]";
% [1] url: https://revistadigital.amca.mx/sdm_downloads/a-note-on-the-energy-shaping-control-applied-to-mechanical-systems-a-pendulum-example/


% Experiment: 
% Based on the lecture where the conditions for ESPACILAG is equivalent to IDA-PBC were taught:
    % \alpha = q
    % \phi = 0
    % q_a= q
    % W = 1
    % V_d = U_a
    % J_a = (M / Md) * (M / Md) * Jd
    % Da = (M / Md) * (M / Md) * Rd


% Task:
% Implementing the ESPACILAG based controller equivalent to the IDA-PBC controller with unbounded torques for
% the Inertia Wheel Pendulum from: 
    % "Control of the Inertia Wheel Pendulum by Bounded Torques" [2], derived from
    % "Stabilization of a class of underactuated mechanical systems via interconnection and damping assignment" [3].  
    
    % [2] doi: https://doi.org/10.1016/j.ifacol.2020.12.1749
    % [3] doi: 10.1109/TAC.2002.800770

function ESPACILAG_IWP_IDA_PBC_qp_nbt()
    % --- Physical Parameters ---
    m3 = 10; 
    I1 = 0.1; 
    I2 = 0.2; 
    u_max = 45; 
    
    % Open-loop Inertia Matrix
    M = [I1, 0; 
         0,  I2]; 
     
    % Derived/Chosen parameters and Controller Gains.
    a1 = 1; 
    a2 = -1.5; 
    a3 = 6;
    kp = 3.75;
    kv = 10;
    
    gamma1 = 30;
    gamma2 = 4.5;
    
    Md = [a1, a2; 
          a2, a3]; 
      
    delta = a1*a3 - a2^2;
    k2 = -(I2 * (a1 + a2)) / delta; 
    k1 = kp * k2; 
    
    % ESPACILAG Transformations
    Ma = M * (Md \ M); 
    Ta = M / Md;          
    Ta_inv = Md / M;      
    
    G = [-1; 1];
    GT = G.';
    Rd_p = G * kv * GT;
    
    Da = Ta * Rd_p * Ta.'; 
    
    % Pack parameters into a structure for the ODE function
    params = struct('M', M, 'I1', I1, 'I2', I2, 'm3', m3, 'G', G, ...
                    'Ta_inv', Ta_inv, 'Da', Da, 'gamma2', gamma2, 'k1', k1);

    % Initial configurations from Section V
    % Rest initial state: [q1(0); q2(0); p1(0); p2(0)] = [3.14; 0; 0; 0]
    x0 = [3.14; 0; 0; 0]; 
    tspan = [0, 40]; % Time horizon: 40 seconds

    % Run simulation    
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, x] = ode45(@(t, x) system_dynamics(t, x, params), tspan, x0, options);

    % Compute the control torque post-simulation for plotting
    u_computed = zeros(length(t), 1);
    for i = 1:length(t)
        [~, u_computed(i)] = compute_espacilag_tau(x(i, :)', params);
    end

    % --- Plotting Results ---
    
    % Figure 3 Replication: Joint Positions q1 and q2
    figure('Position', [100, 100, 550, 400]);
    plot(t, x(:, 1), 'k-', 'LineWidth', 1.5); hold on;
    plot(t, x(:, 2), 'k--', 'LineWidth', 1.2);
    grid on; xlabel('time (seconds)'); ylabel('Joint Positions [rad]');
    legend('q_1', 'q_2'); 
    title('ESPACILAG (q, p): Replication of Figure 3');
    ylim([-2, 3.5]);

    % Figure 5 Replication: Applied Non-Bounded Torque u
    figure('Position', [680, 100, 550, 400]);
    plot(t, u_computed, 'k-', 'LineWidth', 1.5); hold on;
    % Plot the u_max threshold boundary for visualization
    plot(tspan, [u_max, u_max], 'r:', 'LineWidth', 1.2);
    plot(tspan, [-u_max, -u_max], 'r:', 'LineWidth', 1.2);
    grid on; xlabel('time (seconds)'); ylabel('Control Torque (Nm)');
    title('ESPACILAG (q, p): Replication of Figure 5.');
    ylim([-55, 55]);
    
    % Display peak torque value achieved
    fprintf('Peak absolute torque achieved: %.2f Nm\n', max(abs(u_computed)));
end

% Modular Torque Computation Function
function [tau, u] = compute_espacilag_tau(x, params)
    q = x(1:2);
    p = x(3:4);
    
    m3 = params.m3;
    I1 = params.I1;
    I2 = params.I2;
    Ta_inv = params.Ta_inv;
    Da = params.Da;
    gamma2 = params.gamma2;
    k1 = params.k1;
    
    a1 = Ta_inv(1,1)*I1; 
    a2 = Ta_inv(2,1)*I1;

    % 1. Natural Plant Gravitational Vector
    grad_q_H = [-m3 * sin(q(1)); 0]; 

    % Open loop momentum gradient: exactly matches grad_pa_Ha
    grad_p_H = [p(1) / I1; 
                p(2) / I2];

    % 2. Precise Analytical ESPACILAG Potential Energy Gradient Vector
    dVd_dq1 = -(I1 * m3 / (a1 + a2)) * sin(q(1)) + k1 * gamma2 * (q(2) + gamma2 * q(1));
    dVd_dq2 = k1 * (q(2) + gamma2 * q(1));
    grad_qa_Ha = [dVd_dq1; dVd_dq2];

    % 3. Energy Shaping Torque Vector
    tau_es = grad_q_H - Ta_inv * grad_qa_Ha; 

    % 4. Damping Injection Torque Vector (Mapped via matrices)
    tau_di = -Ta_inv * Da * grad_p_H;
    
    % Total generalized control torque vector
    tau = tau_es + tau_di; 
    
    % Project tau into the scalar actuator input space (u)
    u = tau(2); 
end

% Desired closed loop acquisition:
function dxdt = system_dynamics(~, x, params)
    q1 = x(1);
    p1 = x(3);
    p2 = x(4);
    
    I1 = params.I1;
    I2 = params.I2;
    m3 = params.m3;

    % Get control input u from the shared controller function
    [~, u] = compute_espacilag_tau(x, params);

    % Hamiltonian Equations of Motion
    dq1_dt = p1 / I1;
    dq2_dt = p2 / I2;
    dp1_dt = m3 * sin(q1) - u;
    dp2_dt = u;

    dxdt = [dq1_dt; dq2_dt; dp1_dt; dp2_dt];
end