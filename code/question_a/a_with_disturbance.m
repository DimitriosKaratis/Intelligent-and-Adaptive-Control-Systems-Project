%% ΕΥΦΥΗ ΚΑΙ ΠΡΟΣΑΡΜΟΣΤΙΚΑ ΣΥΣΤΗΜΑΤΑ ΑΥΤΟΜΑΤΟΥ ΕΛΕΓΧΟΥ - Dec 2025
% ΚΑΡΑΤΗΣ ΔΗΜΗΤΡΙΟΣ 10775 

%% ==========================================================================
%  Προσομοίωση Ελεγκτή Γραμμικοποίησης (Ερώτημα α) - ΜΕ ΔΙΑΤΑΡΑΧΗ
%  ==========================================================================

clc; clear; close all;

% Ορισμός Παραμέτρων Ελέγχου
theta_star = 2;
c1 = 10;
c2 = 20;

% Ορισμός Παραμέτρων Διαταραχής 
T_FINAL = 6;                  % Συνολικός χρόνος προσομοίωσης
T_DISTURB_START = 3.0;        % Χρόνος έναρξης διαταραχής
DISTURBANCE_MAG = 50;         % Πλατος διαταραχής
DISTURBANCE_DURATION = 0.07;  % Διάρκεια διαταραχής σε seconds

% Ορισμός Αρχικών Συνθηκών
X0_cases = [
    0.5, 0;    % Case 1
    -2, 20;    % Case 2
    2, 20;     % Case 3
    -1, 10     % Case 4
];

num_cases = size(X0_cases, 1);
t_span = [0 T_FINAL];         % Χρόνος προσομοίωσης
LINE_WIDTH = 2;

% Ορισμός της Συνάρτησης Δυναμικής (ODE) ΜΕ ΔΙΑΤΑΡΑΧΗ
f = @(t, X) [
    X(2) + theta_star * X(1)^3; % x1_dot
    (-c2 * X(1) - c1 * X(2) - c1 * theta_star * X(1)^3 - theta_star * 3*X(1)^2 * X(2) - theta_star^2 * 3*X(1)^5) + ... % u
    (t >= T_DISTURB_START && t <= T_DISTURB_START + DISTURBANCE_DURATION) * DISTURBANCE_MAG % + d(t)
];

% Εκτέλεση Προσομοιώσεων
h_x1 = figure('Name', 'Γράφημα 1: Μεταβλητή x1 - ΜΕ ΔΙΑΤΑΡΑΧΗ');
h_x2 = figure('Name', 'Γράφημα 2: Μεταβλητή x2 - ΜΕ ΔΙΑΤΑΡΑΧΗ');
h_u_total = figure('Name', 'Γράφημα 3: Είσοδος Ελέγχου (u + d)'); % Συνολική Είσοδος

figure(h_x1); hold on; grid on; title(['Απόκριση x_1(t) σε Εξωτερική Διαταραχή (' num2str(T_DISTURB_START) 's)']); xlabel('Χρόνος t (s)'); ylabel('x_1');
figure(h_x2); hold on; grid on; title(['Απόκριση x_2(t) σε Εξωτερική Διαταραχή (' num2str(T_DISTURB_START) 's)']); xlabel('Χρόνος t (s)'); ylabel('x_2');
figure(h_u_total); hold on; grid on; title('Συνολική Είσοδος Ελέγχου (u_{total} = u + d)'); xlabel('Χρόνος t (s)'); ylabel('u_{total}');


% Βρόχος Εκτέλεσης
for i = 1:num_cases
    X0 = X0_cases(i, :);
    
    % Επίλυση ODE
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, X] = ode15s(f, t_span, X0, options);
    
    % Υπολογισμός Ελεγκτή u(t) και Διαταραχής d(t)
    x1 = X(:, 1);
    x2 = X(:, 2);
    g_x1 = x1.^3;
    g_prime_x1 = 3*x1.^2;
    
    u = -c2*x1 - c1*x2 - c1*theta_star*g_x1 - theta_star.*g_prime_x1.*x2 - theta_star^2.*g_prime_x1.*g_x1;
    
    % Υπολογισμός Διαταραχής d(t)
    dist = (t >= T_DISTURB_START & t <= T_DISTURB_START + DISTURBANCE_DURATION) * DISTURBANCE_MAG;
    u_total = u + dist;
    
    % Figures
    plot_name_case = ['Case ' num2str(i) ' (x_1(0)=' num2str(X0(1)) ', x_2(0)=' num2str(X0(2)) ')'];

    figure(h_x1);
    plot(t, x1, 'LineWidth', LINE_WIDTH, 'DisplayName', plot_name_case);
    
    figure(h_x2);
    plot(t, x2, 'LineWidth', LINE_WIDTH, 'DisplayName', plot_name_case);
    
    figure(h_u_total);
    plot(t, u_total, 'LineWidth', LINE_WIDTH, 'DisplayName', plot_name_case);
end

% Προσθήκη Legends σε όλα τα Figures
figure(h_x1); legend('show', 'Location', 'northeast'); 
figure(h_x2); legend('show', 'Location', 'southeast'); 
figure(h_u_total); legend('show', 'Location', 'southeast');

% Αποθήκευση Γραφημάτων
saveas(h_x1, 'x1_plot_with_disturbance.png'); 
saveas(h_x2, 'x2_plot_with_disturbance.png'); 
saveas(h_u_total, 'u_plot_with_disturbance.png');  

hold off;
