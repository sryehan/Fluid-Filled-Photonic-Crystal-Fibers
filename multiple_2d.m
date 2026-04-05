%% S1 Scenario - Bright Soliton - Oscilloscope Style Plots
clear all; close all; clc;

% Parameters - CHANGE ω VALUES HERE AS YOU WANT
w_vals = [0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0];  % You can add/remove any values
% w_vals = 0.5:0.3:3;  % Alternative: automatic range from 0.5 to 3 with step 0.3
% w_vals = [0.5, 1, 1.5, 2, 2.5];  % Or specific values

m=0.5; n=-0.8; d=0.1; g=0.2; c=-0.5;
a=1; b=1; r=0.5; e1=1; s1=1; p0=0.5; p1=0.5; i1=0.5;
l=0; rho=0.01;

alpha = m^2 + 2*(n + d*g);
disp(['Alpha = ', num2str(alpha)]);
disp(['Number of ω values = ', num2str(length(w_vals))]);

% Generate colors automatically for as many curves as needed
colors = jet(length(w_vals));  % Using jet colormap for many distinct colors
% Alternative colormaps: hsv, parula, turbo, lines

% Create folder
if ~exist('soliton_figures_oscilloscope', 'dir')
    mkdir('soliton_figures_oscilloscope');
end

%% Spatial profiles (x-axis, t=0)
fprintf('Calculating spatial profiles (t=0) for %d ω values...\n', length(w_vals));

% Space
x = linspace(-8, 8, 1000);

% Preallocate
L1_x = zeros(length(w_vals), length(x));
L4_x = zeros(length(w_vals), length(x));
L5_x = zeros(length(w_vals), length(x));

for k = 1:length(w_vals)
    w = w_vals(k);
    fprintf('  Computing ω = %.3f...\n', w);
    
    for i = 1:length(x)
        xx = x(i);
        tt = 0;  % Fixed time
        
        % L1 calculation
        theta = w*(xx - m*tt);
        sh = sinh(b*theta);
        ch = cosh(b*theta);
        phi = m*xx + tt/2*(b^2*w^2 - m^2 - 2*d*g) + l + rho*randn;
        den = c*s1^2 + 4*a^2*r^2*w^2*e1^2*(1 + 2*ch.*sh - 2*ch.^2);
        if abs(den)>1e-6
            L1_x(k,i) = abs((4*a*b*r*w^2*s1*e1*(sh - ch))/den * exp(1i*phi));
        end
        
        % L4 calculation
        theta_L4 = (s1*sqrt(c))/(a*p0)*(xx - m*tt);
        phi_L4 = m*xx - tt/(4*a^2*p0^2)*(b^2*c*s1^2 + 2*a^2*m^2*p0^2 + 4*a^2*d*g*p0^2) + l + rho*randn;
        num_L4 = b*s1*(b^2*r^2 - a^2 + 2*a*b*r*sinh(b*theta_L4));
        den_L4 = 2*a*p0*(a^2 + b^2*r^2 - 2*a*b*r*cosh(b*theta_L4));
        if abs(den_L4)>1e-6
            L4_x(k,i) = abs(num_L4/den_L4 * exp(1i*phi_L4));
        end
        
        % L5 calculation
        theta_L5 = (2*i1*sqrt(c))/(b*e1)*(xx - m*tt);
        phi_L5 = m*xx - tt/(2*e1^2)*(m^2*e1^2 + 2*c*i1^2 + 2*d*g*e1^2) + l + rho*randn;
        term1 = sinh(b*theta_L5) - cosh(b*theta_L5);
        term2 = cosh(b*theta_L5) - sinh(b*theta_L5);
        num_L5 = i1*(a*(i1*p1 - s1*e1) + b*r*(i1*p1 + s1*e1)*term1);
        den_L5 = e1*(a*(i1*p1 - s1*e1) + b*r*(i1*p1 + s1*e1)*term2);
        if abs(den_L5)>1e-6
            L5_x(k,i) = abs(-num_L5/den_L5 * exp(1i*phi_L5));
        end
    end
    
    % Normalize each profile
    L1_x(k,:) = L1_x(k,:) / max(L1_x(k,:));
    L4_x(k,:) = L4_x(k,:) / max(L4_x(k,:));
    L5_x(k,:) = L5_x(k,:) / max(L5_x(k,:));
end

%% Temporal profiles (t-axis, x=0)
fprintf('Calculating temporal profiles (x=0) for %d ω values...\n', length(w_vals));

% Time
t = linspace(-3, 3, 1000);

% Preallocate
L1_t = zeros(length(w_vals), length(t));
L4_t = zeros(length(w_vals), length(t));
L5_t = zeros(length(w_vals), length(t));

for k = 1:length(w_vals)
    w = w_vals(k);
    fprintf('  Computing ω = %.3f...\n', w);
    
    for j = 1:length(t)
        xx = 0;  % Fixed position
        tt = t(j);
        
        % L1 calculation
        theta = w*(xx - m*tt);
        sh = sinh(b*theta);
        ch = cosh(b*theta);
        phi = m*xx + tt/2*(b^2*w^2 - m^2 - 2*d*g) + l + rho*randn;
        den = c*s1^2 + 4*a^2*r^2*w^2*e1^2*(1 + 2*ch.*sh - 2*ch.^2);
        if abs(den)>1e-6
            L1_t(k,j) = abs((4*a*b*r*w^2*s1*e1*(sh - ch))/den * exp(1i*phi));
        end
        
        % L4 calculation
        theta_L4 = (s1*sqrt(c))/(a*p0)*(xx - m*tt);
        phi_L4 = m*xx - tt/(4*a^2*p0^2)*(b^2*c*s1^2 + 2*a^2*m^2*p0^2 + 4*a^2*d*g*p0^2) + l + rho*randn;
        num_L4 = b*s1*(b^2*r^2 - a^2 + 2*a*b*r*sinh(b*theta_L4));
        den_L4 = 2*a*p0*(a^2 + b^2*r^2 - 2*a*b*r*cosh(b*theta_L4));
        if abs(den_L4)>1e-6
            L4_t(k,j) = abs(num_L4/den_L4 * exp(1i*phi_L4));
        end
        
        % L5 calculation
        theta_L5 = (2*i1*sqrt(c))/(b*e1)*(xx - m*tt);
        phi_L5 = m*xx - tt/(2*e1^2)*(m^2*e1^2 + 2*c*i1^2 + 2*d*g*e1^2) + l + rho*randn;
        term1 = sinh(b*theta_L5) - cosh(b*theta_L5);
        term2 = cosh(b*theta_L5) - sinh(b*theta_L5);
        num_L5 = i1*(a*(i1*p1 - s1*e1) + b*r*(i1*p1 + s1*e1)*term1);
        den_L5 = e1*(a*(i1*p1 - s1*e1) + b*r*(i1*p1 + s1*e1)*term2);
        if abs(den_L5)>1e-6
            L5_t(k,j) = abs(-num_L5/den_L5 * exp(1i*phi_L5));
        end
    end
    
    % Normalize each profile
    L1_t(k,:) = L1_t(k,:) / max(L1_t(k,:));
    L4_t(k,:) = L4_t(k,:) / max(L4_t(k,:));
    L5_t(k,:) = L5_t(k,:) / max(L5_t(k,:));
end

%% FIGURE 1: Spatial Profiles for L1 (t=0)
fig1 = figure('Position', [100, 100, 1100, 700], 'Color', 'k');
hold on;
for k = 1:length(w_vals)
    plot(x, L1_x(k,:), 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('ω = %.3f', w_vals(k)));
end
hold off;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.3, 0.3, 0.3]);
grid on;
xlabel('Spatial Coordinate (x) →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Intensity |L₁| →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('L₁ Equation - Spatial Profile (t = 0) [%d ω values]', length(w_vals)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'TextColor', 'w', 'Color', 'k', 'EdgeColor', 'w', 'FontSize', 8);
xlim([-8, 8]);
ylim([0, 1.1]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
saveas(fig1, 'soliton_figures_oscilloscope/Figure1_L1_Spatial_Profile.png');
saveas(fig1, 'soliton_figures_oscilloscope/Figure1_L1_Spatial_Profile.fig');
fprintf('✓ Figure 1 saved: L1_Spatial_Profile (%d curves)\n', length(w_vals));

%% FIGURE 2: Spatial Profiles for L4 (t=0)
fig2 = figure('Position', [150, 150, 1100, 700], 'Color', 'k');
hold on;
for k = 1:length(w_vals)
    plot(x, L4_x(k,:), 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('ω = %.3f', w_vals(k)));
end
hold off;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.3, 0.3, 0.3]);
grid on;
xlabel('Spatial Coordinate (x) →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Intensity |L₄| →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('L₄ Equation - Spatial Profile (t = 0) [%d ω values]', length(w_vals)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'TextColor', 'w', 'Color', 'k', 'EdgeColor', 'w', 'FontSize', 8);
xlim([-8, 8]);
ylim([0, 1.1]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
saveas(fig2, 'soliton_figures_oscilloscope/Figure2_L4_Spatial_Profile.png');
saveas(fig2, 'soliton_figures_oscilloscope/Figure2_L4_Spatial_Profile.fig');
fprintf('✓ Figure 2 saved: L4_Spatial_Profile (%d curves)\n', length(w_vals));

%% FIGURE 3: Spatial Profiles for L5 (t=0)
fig3 = figure('Position', [200, 200, 1100, 700], 'Color', 'k');
hold on;
for k = 1:length(w_vals)
    plot(x, L5_x(k,:), 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('ω = %.3f', w_vals(k)));
end
hold off;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.3, 0.3, 0.3]);
grid on;
xlabel('Spatial Coordinate (x) →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Intensity |L₅| →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('L₅ Equation - Spatial Profile (t = 0) [%d ω values]', length(w_vals)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'TextColor', 'w', 'Color', 'k', 'EdgeColor', 'w', 'FontSize', 8);
xlim([-8, 8]);
ylim([0, 1.1]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
saveas(fig3, 'soliton_figures_oscilloscope/Figure3_L5_Spatial_Profile.png');
saveas(fig3, 'soliton_figures_oscilloscope/Figure3_L5_Spatial_Profile.fig');
fprintf('✓ Figure 3 saved: L5_Spatial_Profile (%d curves)\n', length(w_vals));

%% FIGURE 4: Temporal Profiles for L1 (x=0)
fig4 = figure('Position', [100, 100, 1100, 700], 'Color', 'k');
hold on;
for k = 1:length(w_vals)
    plot(t, L1_t(k,:), 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('ω = %.3f', w_vals(k)));
end
hold off;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.3, 0.3, 0.3]);
grid on;
xlabel('Time (t) →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Intensity |L₁| →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('L₁ Equation - Temporal Profile (x = 0) [%d ω values]', length(w_vals)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'TextColor', 'w', 'Color', 'k', 'EdgeColor', 'w', 'FontSize', 8);
xlim([-3, 3]);
ylim([0, 1.1]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
saveas(fig4, 'soliton_figures_oscilloscope/Figure4_L1_Temporal_Profile.png');
saveas(fig4, 'soliton_figures_oscilloscope/Figure4_L1_Temporal_Profile.fig');
fprintf('✓ Figure 4 saved: L1_Temporal_Profile (%d curves)\n', length(w_vals));

%% FIGURE 5: Temporal Profiles for L4 (x=0)
fig5 = figure('Position', [150, 150, 1100, 700], 'Color', 'k');
hold on;
for k = 1:length(w_vals)
    plot(t, L4_t(k,:), 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('ω = %.3f', w_vals(k)));
end
hold off;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.3, 0.3, 0.3]);
grid on;
xlabel('Time (t) →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Intensity |L₄| →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('L₄ Equation - Temporal Profile (x = 0) [%d ω values]', length(w_vals)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'TextColor', 'w', 'Color', 'k', 'EdgeColor', 'w', 'FontSize', 8);
xlim([-3, 3]);
ylim([0, 1.1]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
saveas(fig5, 'soliton_figures_oscilloscope/Figure5_L4_Temporal_Profile.png');
saveas(fig5, 'soliton_figures_oscilloscope/Figure5_L4_Temporal_Profile.fig');
fprintf('✓ Figure 5 saved: L4_Temporal_Profile (%d curves)\n', length(w_vals));

%% FIGURE 6: Temporal Profiles for L5 (x=0)
fig6 = figure('Position', [200, 200, 1100, 700], 'Color', 'k');
hold on;
for k = 1:length(w_vals)
    plot(t, L5_t(k,:), 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', sprintf('ω = %.3f', w_vals(k)));
end
hold off;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', [0.3, 0.3, 0.3]);
grid on;
xlabel('Time (t) →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Intensity |L₅| →', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('L₅ Equation - Temporal Profile (x = 0) [%d ω values]', length(w_vals)), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'TextColor', 'w', 'Color', 'k', 'EdgeColor', 'w', 'FontSize', 8);
xlim([-3, 3]);
ylim([0, 1.1]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'Box', 'on');
saveas(fig6, 'soliton_figures_oscilloscope/Figure6_L5_Temporal_Profile.png');
saveas(fig6, 'soliton_figures_oscilloscope/Figure6_L5_Temporal_Profile.fig');
fprintf('✓ Figure 6 saved: L5_Temporal_Profile (%d curves)\n', length(w_vals));

fprintf('\n========================================\n');
fprintf('✓ ALL 6 OSCILLOSCOPE-STYLE FIGURES SAVED\n');
fprintf('  in "soliton_figures_oscilloscope" folder\n');
fprintf('========================================\n');
fprintf('Total ω values used: %d\n', length(w_vals));
fprintf('ω range: %.3f to %.3f\n', min(w_vals), max(w_vals));
fprintf('\nFIGURES (each has %d curves):\n', length(w_vals));
fprintf('  1. L1 - Spatial Profile (t=0, x-axis)\n');
fprintf('  2. L4 - Spatial Profile (t=0, x-axis)\n');
fprintf('  3. L5 - Spatial Profile (t=0, x-axis)\n');
fprintf('  4. L1 - Temporal Profile (x=0, t-axis)\n');
fprintf('  5. L4 - Temporal Profile (x=0, t-axis)\n');
fprintf('  6. L5 - Temporal Profile (x=0, t-axis)\n');
fprintf('========================================\n');
disp('✨ Complete! ✨');