%% ΕΥΦΥΗ ΚΑΙ ΠΡΟΣΑΡΜΟΣΤΙΚΑ ΣΥΣΤΗΜΑΤΑ ΑΥΤΟΜΑΤΟΥ ΕΛΕΓΧΟΥ - Dec 2025
% ΚΑΡΑΤΗΣ ΔΗΜΗΤΡΙΟΣ 10775

%% ==========================================================================
%  Προσομοίωση Ελεγκτή Γραμμικοποίησης (Ερώτημα α)
%  ==========================================================================

clc; clear; close all;

% Ορισμός Παραμέτρων
theta_star = 2;
c1 = 10; 
c2 = 20; 

% Ορισμός Αρχικών Συνθηκών
X0_cases = [
    0.5, 0;        % Case 1
    -2, 20;        % Case 2
    2, 20;         % Case 3
    -1, 10         % Case 4
];

num_cases = size(X0_cases, 1);
t_span = [0 2.5];  % Χρόνος προσομοίωσης
LINE_WIDTH = 2;  

% Ορισμός της Συνάρτησης Δυναμικής (ODE)
f = @(t, X) [
    X(2) + theta_star * X(1)^3;     % x1_dot = x2 + theta*g(x1)
    -c2 * X(1) - c1 * X(2) - c1 * theta_star * X(1)^3 - theta_star * 3*X(1)^2 * X(2) ...
    - theta_star^2 * 3*X(1)^5       % x2_dot = u
];

% Εκτέλεση Προσομοιώσεων 
h_x1 = figure('Name', 'Γράφημα 1: Μεταβλητή x1');
h_x2 = figure('Name', 'Γράφημα 2: Μεταβλητή x2');
h_u = figure('Name', 'Γράφημα 3: Είσοδος Ελέγχου u');

figure(h_x1); hold on; grid on; title('Μεταβλητή Κατάστασης x_1(t)'); xlabel('Χρόνος t (s)'); ylabel('x_1');
figure(h_x2); hold on; grid on; title('Μεταβλητή Κατάστασης x_2(t)'); xlabel('Χρόνος t (s)'); ylabel('x_2');
figure(h_u); hold on; grid on; title('Είσοδος Ελέγχου u(t)'); xlabel('Χρόνος t (s)'); ylabel('u');

% Βρόχος Εκτέλεσης
for i = 1:num_cases
    X0 = X0_cases(i, :);
    
    % Επίλυση ODE
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, X] = ode15s(f, t_span, X0);
    
    % Υπολογισμός Ελεγκτή u(t)
    x1 = X(:, 1);
    x2 = X(:, 2);
    g_x1 = x1.^3;
    g_prime_x1 = 3*x1.^2;
    
    u = -c2*x1 - c1*x2 - c1*theta_star*g_x1 - theta_star.*g_prime_x1.*x2 - theta_star^2.*g_prime_x1.*g_x1;
        
    % Γράφημα 1: x1
    figure(h_x1);
    plot(t, x1, 'LineWidth', LINE_WIDTH, 'DisplayName', ['Case ' num2str(i) ' (x_1(0)=' num2str(X0(1)) ', x_2(0)=' num2str(X0(2)) ')']);
    
    % Γράφημα 2: x2
    figure(h_x2);
    plot(t, x2, 'LineWidth', LINE_WIDTH, 'DisplayName', ['Case ' num2str(i) ' (x_1(0)=' num2str(X0(1)) ', x_2(0)=' num2str(X0(2)) ')']);
    
    % Γράφημα 3: u
    figure(h_u);
    plot(t, u, 'LineWidth', LINE_WIDTH, 'DisplayName', ['Case ' num2str(i) ' (x_1(0)=' num2str(X0(1)) ', x_2(0)=' num2str(X0(2)) ')']);
end


% Προσθήκη Legends σε όλα τα Figures
figure(h_x1); legend('show', 'Location', 'northeast');
figure(h_x2); legend('show', 'Location', 'southeast'); 
figure(h_u); legend('show', 'Location', 'southeast');  

% Αποθήκευση Γραφημάτων
saveas(h_x1, 'x1_plot.png');
saveas(h_x2, 'x2_plot.png'); 
saveas(h_u, 'u_plot.png'); 

hold off;
