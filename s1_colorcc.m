%% S1 Scenario - Bright Soliton - Custom Hex Color Palette
clear all; close all; clc;

% Parameters
w=1; m=0.5; n=-0.8; d=0.1; g=0.2; c=-0.5;
a=1; b=1; r=0.5; e1=1; s1=1; p0=0.5; p1=0.5; i1=0.5;
l=0; rho=0.01;

alpha = m^2 + 2*(n + d*g);
disp(['Alpha = ', num2str(alpha)]);

% Space and time
x = linspace(-8,8,500);
t = linspace(-3,3,400);
[X,T] = meshgrid(x,t);

% Preallocate for L1, L4, L5
L1 = zeros(size(X));
L4 = zeros(size(X));
L5 = zeros(size(X));

% Calculate L1, L4, L5
fprintf('Calculating...\n');
for i=1:length(x)
    for j=1:length(t)
        xx = X(j,i);
        tt = T(j,i);
        
        % L1 calculation
        theta = w*(xx - m*tt);
        sh = sinh(b*theta);
        ch = cosh(b*theta);
        
        phi = m*xx + tt/2*(b^2*w^2 - m^2 - 2*d*g) + l + rho*randn;
        den = c*s1^2 + 4*a^2*r^2*w^2*e1^2*(1 + 2*ch.*sh - 2*ch.^2);
        if abs(den)>1e-6
            L1(j,i) = abs((4*a*b*r*w^2*s1*e1*(sh - ch))/den * exp(1i*phi));
        end
        
        % L4 calculation
        theta_L4 = (s1*sqrt(c))/(a*p0)*(xx - m*tt);
        phi_L4 = m*xx - tt/(4*a^2*p0^2)*(b^2*c*s1^2 + 2*a^2*m^2*p0^2 + 4*a^2*d*g*p0^2) + l + rho*randn;
        num_L4 = b*s1*(b^2*r^2 - a^2 + 2*a*b*r*sinh(b*theta_L4));
        den_L4 = 2*a*p0*(a^2 + b^2*r^2 - 2*a*b*r*cosh(b*theta_L4));
        if abs(den_L4)>1e-6
            L4(j,i) = abs(num_L4/den_L4 * exp(1i*phi_L4));
        end
        
        % L5 calculation
        theta_L5 = (2*i1*sqrt(c))/(b*e1)*(xx - m*tt);
        phi_L5 = m*xx - tt/(2*e1^2)*(m^2*e1^2 + 2*c*i1^2 + 2*d*g*e1^2) + l + rho*randn;
        term1 = sinh(b*theta_L5) - cosh(b*theta_L5);
        term2 = cosh(b*theta_L5) - sinh(b*theta_L5);
        num_L5 = i1*(a*(i1*p1 - s1*e1) + b*r*(i1*p1 + s1*e1)*term1);
        den_L5 = e1*(a*(i1*p1 - s1*e1) + b*r*(i1*p1 + s1*e1)*term2);
        if abs(den_L5)>1e-6
            L5(j,i) = abs(-num_L5/den_L5 * exp(1i*phi_L5));
        end
    end
end

% Normalize
L1_norm = L1 / max(L1(:));
L4_norm = L4 / max(L4(:));
L5_norm = L5 / max(L5(:));

% Create custom colormap from your hex colors
hex_colors = {'#473335', '#548687', '#ffffc7', '#b0413e', '#fcaa67'};

% Convert hex to RGB
rgb_colors = zeros(length(hex_colors), 3);
for i = 1:length(hex_colors)
    hex = hex_colors{i};
    rgb_colors(i,1) = hex2dec(hex(2:3))/255;
    rgb_colors(i,2) = hex2dec(hex(4:5))/255;
    rgb_colors(i,3) = hex2dec(hex(6:7))/255;
end

% Create interpolated colormap
n_interp = 256;
custom_cmap = zeros(n_interp, 3);
for i = 1:3
    custom_cmap(:,i) = interp1(linspace(0,1,length(rgb_colors)), rgb_colors(:,i), linspace(0,1,n_interp));
end

% Create folder
if ~exist('soliton_figures', 'dir')
    mkdir('soliton_figures');
end

%% FIGURE 1: L1 - 3D Surface (Separate Figure)
fig1 = figure('Position', [100, 100, 900, 700], 'Color', rgb_colors(1,:)*0.3);
surf(X, T, L1_norm, 'EdgeColor', 'none', 'FaceAlpha', 0.95);
colormap(custom_cmap);
view(50, 35);
title('L₁ (Eq 1) - Bright Soliton Propagation', 'Color', hex2rgb('#ffffc7'), 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spatial Coordinate (x) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
ylabel('Time (t) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
zlabel('Intensity |L₁| →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
set(gca, 'Color', rgb_colors(1,:)*0.3, 'XColor', hex2rgb('#ffffc7'), 'YColor', hex2rgb('#ffffc7'), 'ZColor', hex2rgb('#ffffc7'));
grid on;
set(gca, 'GridColor', hex2rgb('#548687'), 'GridAlpha', 0.3);
material([0.6, 0.8, 0.4]);
lighting gouraud;
light('Position', [2, 1, 3], 'Color', hex2rgb('#fcaa67'));
light('Position', [-1, -1, 2], 'Color', hex2rgb('#b0413e'));
colorbar('Color', hex2rgb('#ffffc7'));
saveas(fig1, 'soliton_figures/Figure1_L1_3D.png');
saveas(fig1, 'soliton_figures/Figure1_L1_3D.fig');
fprintf('✓ Figure 1 saved: L1_3D\n');

%% FIGURE 2: L4 - 3D Surface (Separate Figure)
fig2 = figure('Position', [150, 150, 900, 700], 'Color', rgb_colors(1,:)*0.3);
surf(X, T, L4_norm, 'EdgeColor', 'none', 'FaceAlpha', 0.95);
colormap(custom_cmap);
view(50, 35);
title('L₄ (Eq 4) - Bright Soliton Propagation', 'Color', hex2rgb('#ffffc7'), 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spatial Coordinate (x) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
ylabel('Time (t) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
zlabel('Intensity |L₄| →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
set(gca, 'Color', rgb_colors(1,:)*0.3, 'XColor', hex2rgb('#ffffc7'), 'YColor', hex2rgb('#ffffc7'), 'ZColor', hex2rgb('#ffffc7'));
grid on;
set(gca, 'GridColor', hex2rgb('#548687'), 'GridAlpha', 0.3);
material([0.6, 0.8, 0.4]);
lighting gouraud;
light('Position', [2, 1, 3], 'Color', hex2rgb('#fcaa67'));
light('Position', [-1, -1, 2], 'Color', hex2rgb('#b0413e'));
colorbar('Color', hex2rgb('#ffffc7'));
saveas(fig2, 'soliton_figures/Figure2_L4_3D.png');
saveas(fig2, 'soliton_figures/Figure2_L4_3D.fig');
fprintf('✓ Figure 2 saved: L4_3D\n');

%% FIGURE 3: L5 - 3D Surface (Separate Figure)
fig3 = figure('Position', [200, 200, 900, 700], 'Color', rgb_colors(1,:)*0.3);
surf(X, T, L5_norm, 'EdgeColor', 'none', 'FaceAlpha', 0.95);
colormap(custom_cmap);
view(50, 35);
title('L₅ (Eq 5) - Bright Soliton Propagation', 'Color', hex2rgb('#ffffc7'), 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spatial Coordinate (x) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
ylabel('Time (t) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
zlabel('Intensity |L₅| →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
set(gca, 'Color', rgb_colors(1,:)*0.3, 'XColor', hex2rgb('#ffffc7'), 'YColor', hex2rgb('#ffffc7'), 'ZColor', hex2rgb('#ffffc7'));
grid on;
set(gca, 'GridColor', hex2rgb('#548687'), 'GridAlpha', 0.3);
material([0.6, 0.8, 0.4]);
lighting gouraud;
light('Position', [2, 1, 3], 'Color', hex2rgb('#fcaa67'));
light('Position', [-1, -1, 2], 'Color', hex2rgb('#b0413e'));
colorbar('Color', hex2rgb('#ffffc7'));
saveas(fig3, 'soliton_figures/Figure3_L5_3D.png');
saveas(fig3, 'soliton_figures/Figure3_L5_3D.fig');
fprintf('✓ Figure 3 saved: L5_3D\n');

%% FIGURE 4: L1 - 2D Intensity Map (Separate Figure)
fig4 = figure('Position', [100, 100, 800, 600], 'Color', rgb_colors(1,:)*0.3);
imagesc(x, t, L1_norm);
colormap(custom_cmap);
set(gca, 'YDir', 'normal', 'Color', rgb_colors(1,:)*0.3, 'XColor', hex2rgb('#ffffc7'), 'YColor', hex2rgb('#ffffc7'));
title('L₁ (Eq 1) - Intensity Map', 'Color', hex2rgb('#ffffc7'), 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spatial Coordinate (x) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
ylabel('Time (t) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
colorbar('Color', hex2rgb('#ffffc7'));
set(gca, 'FontSize', 10);
saveas(fig4, 'soliton_figures/Figure4_L1_2D_Map.png');
saveas(fig4, 'soliton_figures/Figure4_L1_2D_Map.fig');
fprintf('✓ Figure 4 saved: L1_2D_Map\n');

%% FIGURE 5: L4 - 2D Intensity Map (Separate Figure)
fig5 = figure('Position', [150, 150, 800, 600], 'Color', rgb_colors(1,:)*0.3);
imagesc(x, t, L4_norm);
colormap(custom_cmap);
set(gca, 'YDir', 'normal', 'Color', rgb_colors(1,:)*0.3, 'XColor', hex2rgb('#ffffc7'), 'YColor', hex2rgb('#ffffc7'));
title('L₄ (Eq 4) - Intensity Map', 'Color', hex2rgb('#ffffc7'), 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spatial Coordinate (x) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
ylabel('Time (t) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
colorbar('Color', hex2rgb('#ffffc7'));
set(gca, 'FontSize', 10);
saveas(fig5, 'soliton_figures/Figure5_L4_2D_Map.png');
saveas(fig5, 'soliton_figures/Figure5_L4_2D_Map.fig');
fprintf('✓ Figure 5 saved: L4_2D_Map\n');

%% FIGURE 6: L5 - 2D Intensity Map (Separate Figure)
fig6 = figure('Position', [200, 200, 800, 600], 'Color', rgb_colors(1,:)*0.3);
imagesc(x, t, L5_norm);
colormap(custom_cmap);
set(gca, 'YDir', 'normal', 'Color', rgb_colors(1,:)*0.3, 'XColor', hex2rgb('#ffffc7'), 'YColor', hex2rgb('#ffffc7'));
title('L₅ (Eq 5) - Intensity Map', 'Color', hex2rgb('#ffffc7'), 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Spatial Coordinate (x) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
ylabel('Time (t) →', 'Color', hex2rgb('#ffffc7'), 'FontSize', 11);
colorbar('Color', hex2rgb('#ffffc7'));
set(gca, 'FontSize', 10);
saveas(fig6, 'soliton_figures/Figure6_L5_2D_Map.png');
saveas(fig6, 'soliton_figures/Figure6_L5_2D_Map.fig');
fprintf('✓ Figure 6 saved: L5_2D_Map\n');

% Helper function to convert hex to RGB
function rgb = hex2rgb(hex)
    hex = char(hex);
    if hex(1) == '#'
        hex = hex(2:end);
    end
    rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))]/255;
end

fprintf('\n========================================\n');
fprintf('✓ ALL 6 FIGURES SAVED in "soliton_figures" folder\n');
fprintf('✓ Figure 1: L1 - 3D Surface\n');
fprintf('✓ Figure 2: L4 - 3D Surface\n');
fprintf('✓ Figure 3: L5 - 3D Surface\n');
fprintf('✓ Figure 4: L1 - 2D Intensity Map\n');
fprintf('✓ Figure 5: L4 - 2D Intensity Map\n');
fprintf('✓ Figure 6: L5 - 2D Intensity Map\n');
fprintf('========================================\n');
disp('✨ Complete! 6 separate figure files created! ✨');