%% ΕΥΦΥΗ ΚΑΙ ΠΡΟΣΑΡΜΟΣΤΙΚΑ ΣΥΣΤΗΜΑΤΑ ΑΥΤΟΜΑΤΟΥ ΕΛΕΓΧΟΥ - Dec 2025
% ΚΑΡΑΤΗΣ ΔΗΜΗΤΡΙΟΣ 10775

%% ================================================================================
%  Προσομοίωση Προσαρμοστικού Backstepping Ελεγκτή ΜΕ ΔΙΑΤΑΡΑΧΗ (Ερώτημα ε)
%  ================================================================================

clc; clear; close all;

% --- Παράμετροι Προσομοίωσης --- 
theta_star = 2;   % Πραγματική (άγνωστη στον ελεγκτή) παράμετρος
gamma = 1;        % Ρυθμός προσαρμογής (Adaptation gain)
k1 = 2;           % Κέρδος ελέγχου βήματος 1
k2 = 2;           % Κέρδος ελέγχου βήματος 2

% --- Παράμετροι Διαταραχής ---
T_DISTURB_START = 5.0;       % Χρόνος έναρξης διαταραχής
DISTURBANCE_MAG = 50;         % Πλάτος διαταραχής
DISTURBANCE_DURATION = 0.07;  % Διάρκεια διαταραχής

% Συνάρτηση g(x1) και η παράγωγός της g'(x1)
g = @(x1) x1^3;   
gp = @(x1) 3*x1^2;

% --- Αρχικές Συνθήκες [x1, x2, theta_hat]  --- 
X0_matrix = [ 0.5,  0,   0;   % Case 1 
             -2.0, 20,   0;   % Case 2 
              2.0, 20,   0;   % Case 3 
             -1.0, 10,   0];  % Case 4 


t_span = [0 8]; % Χρονικό διάστημα προσομοίωσης
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.01);

% --- Loop Προσομοίωσης ---
for j = 1:4
    % Δημιουργία Figure 1 για Cases 1-2 και Figure 2 για Cases 3-4
    if j == 1 || j == 3
        fig_h = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.5]);
        row_idx = 1; 
    else
        row_idx = 2; 
    end
    
    x0 = X0_matrix(j, :);
    
    % Επίλυση με ode15s (περιλαμβάνει χρόνο t για τη διαταραχή)
    [t, states] = ode15s(@(t, x) backstepping_sys_dist(t, x, theta_star, gamma, k1, k2, g, gp, T_DISTURB_START, DISTURBANCE_MAG, DISTURBANCE_DURATION), t_span, x0, options);
    
    x1 = states(:,1); 
    x2 = states(:,2); 
    theta_hat = states(:,3);
    
    u_hist = zeros(length(t), 1);
    for i = 1:length(t)
        z1 = x1(i);
        phi = -k1*z1 - theta_hat(i)*g(z1);
        z2 = x2(i) - phi;
        
        % Νόμος Προσαρμογής
        dot_th = gamma * (z1*g(z1) + z2*(k1 + theta_hat(i)*gp(z1))*g(z1));
        
        % Υπολογισμός διαταραχής d(t) για το γράφημα της συνολικής εισόδου
        dist_t = (t(i) >= T_DISTURB_START && t(i) <= T_DISTURB_START + DISTURBANCE_DURATION) * DISTURBANCE_MAG;
        
        % Νόμος Ελέγχου u (συνολική είσοδος u + d)
        u_hist(i) = -z1 - k2*z2 - (k1 + theta_hat(i)*gp(z1))*(x2(i) + theta_hat(i)*g(z1)) - g(z1)*dot_th + dist_t;
    end


    % --- Subplots ---
    
    % Στήλη 1: States x1, x2
    subplot(2, 3, (row_idx-1)*3 + 1);
    plot(t, x1, 'b', 'LineWidth', 1.5); hold on;
    plot(t, x2, 'r', 'LineWidth', 1.5); 
    xline(T_DISTURB_START, '--k', 'HandleVisibility', 'off'); % Διακεκομμένη γραμμή διαταραχής
    grid on;
    title(['Case ', num2str(j), ': States']);
    ylabel(['x_0=[', num2str(x0(1)), ',', num2str(x0(2)), ']']);
    if row_idx == 1, legend('x_1','x_2'); end
    
    % Στήλη 2: theta_hat
    subplot(2, 3, (row_idx-1)*3 + 2);
    plot(t, theta_hat, 'g', 'LineWidth', 1.5); hold on;
    plot(t, ones(size(t))*theta_star, 'k--', 'LineWidth', 1); 
    xline(T_DISTURB_START, '--k', 'HandleVisibility', 'off'); % Διακεκομμένη γραμμή διαταραχής
    grid on;
    title(['Case ', num2str(j), ': Parameter Estimation']);
    if row_idx == 1, legend('\theta_{hat}','\theta^*'); end
    
    % Στήλη 3: Control Input u (Συνολική είσοδος)
    subplot(2, 3, (row_idx-1)*3 + 3);
    plot(t, u_hist, 'm', 'LineWidth', 1.5); hold on;
    xline(T_DISTURB_START, '--k', 'HandleVisibility', 'off'); % Διακεκομμένη γραμμή διαταραχής
    grid on;
    title(['Case ', num2str(j), ': Total Input (u+d)']);
    if row_idx == 2, xlabel('Time (s)'); end
    
    % --- Αποθήκευση Figures ---
    if j == 2
        saveas(fig_h, 'E_Dist_Cases_1_2.png');
    elseif j == 4
        saveas(fig_h, 'E_Dist_Cases_3_4.png');
    end
end

% Συνάρτηση Δυναμικής Backstepping με Διαταραχή
function dxdt = backstepping_sys_dist(t, x, theta_star, gamma, k1, k2, g, gp, T_S, D_M, D_D)
    x1 = x(1); 
    x2 = x(2); 
    theta_hat = x(3);
    
    % Ορισμός σφαλμάτων z 
    z1 = x1;
    phi = -k1*z1 - theta_hat*g(z1);
    z2 = x2 - phi;
    
    % 1. Νόμος Προσαρμογής 
    dot_theta_hat = gamma * (z1*g(z1) + z2*(k1 + theta_hat*gp(z1))*g(z1));
    
    % 2. Νόμος Ελέγχου u 
    u = -z1 - k2*z2 - (k1 + theta_hat*gp(z1))*(x2 + theta_hat*g(z1)) - g(z1)*dot_theta_hat;
    
    % Εξωτερική Διαταραχή d(t) 
    dist = (t >= T_S && t <= T_S + D_D) * D_M;
    
    % Πραγματική Δυναμική Συστήματος 
    dx1 = x2 + theta_star * g(x1);
    dx2 = u + dist; % Η διαταραχή εισέρχεται στην dot_x2
    
    dxdt = [dx1; dx2; dot_theta_hat];
end