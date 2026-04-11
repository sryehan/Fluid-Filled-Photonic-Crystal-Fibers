% =========================================================================
% PCF Glucose Sensor: Sensitivity & Soliton Analysis
% Author: Shahariar Ryehan
% Journal: Elsevier Q1 (Data-Dependent Soliton Dynamics)
% =========================================================================

clear all; clc; close all;

%% ── GLOBAL STYLE DEFAULTS ───────────────────────────────────────────────
set(0, 'DefaultAxesFontName',  'Times New Roman');
set(0, 'DefaultTextFontName',  'Times New Roman');
set(0, 'DefaultAxesFontSize',  12);
set(0, 'DefaultAxesLineWidth', 1.0);
set(0, 'DefaultLineLineWidth', 1.8);

%% 1. DATA LOADING & PREPROCESSING ────────────────────────────────────────
data       = readmatrix('main_18.txt');
n_core_raw = data(:, 1);
n_eff_raw  = data(:, 2);

[unique_ncore, ~, idx] = unique(n_core_raw);
clean_neff = accumarray(idx, n_eff_raw, [], @max);

%% 2. SENSITIVITY ANALYSIS (Linear Fitting) ───────────────────────────────
p           = polyfit(unique_ncore, clean_neff, 1);
y_fit       = polyval(p, unique_ncore);
Sensitivity = p(1);
R_sq = 1 - sum((clean_neff - y_fit).^2) / sum((clean_neff - mean(clean_neff)).^2);

%% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 1 — SENSITIVITY PLOT
%  Clean Elsevier style: white bg, black axes, red data, blue fit
% ═══════════════════════════════════════════════════════════════════════
fig1 = figure('Units','centimeters','Position',[2 2 16 11],'Color','w');

plot(unique_ncore, clean_neff, 'o', ...
    'Color',[0.80 0.10 0.10], ...
    'MarkerFaceColor',[0.80 0.10 0.10], ...
    'MarkerEdgeColor','w', ...
    'MarkerSize', 7, 'LineWidth', 1.2);
hold on;
plot(unique_ncore, y_fit, '-', ...
    'Color',[0.08 0.35 0.72], 'LineWidth', 2.2);

grid on;
box on;
ax = gca;
ax.GridLineStyle  = ':';
ax.GridAlpha      = 0.35;
ax.GridColor      = [0.5 0.5 0.5];
ax.TickDir        = 'out';
ax.TickLength     = [0.012 0.025];
ax.XMinorTick     = 'on';
ax.YMinorTick     = 'on';

xlabel('Glucose Core Refractive Index, \itn\rm_{core}', 'FontSize',13);
ylabel('Effective Mode Index, \itn\rm_{eff}',           'FontSize',13);
title('Sensitivity Analysis of PCF-Based Glucose Sensor', ...
      'FontSize',14, 'FontWeight','bold');

lg = legend({'COMSOL Simulation Data', ...
    sprintf('Linear Fit: \\itn\\rm_{eff} = %.4f\\cdot\\itn\\rm_{core} + %.4f  (\\itR\\rm^{2} = %.4f)', ...
    p(1), p(2), R_sq)}, ...
    'Location','northwest','FontSize',11,'Box','off');

% Inset annotation box
dim = [0.62 0.20 0.30 0.10];
str = sprintf('\\itS\\rm = %.2e RIU^{-1}', Sensitivity);
annotation('textbox', dim, 'String', str, ...
    'FontSize',12, 'FontWeight','bold', ...
    'BackgroundColor','w', 'EdgeColor',[0.3 0.3 0.3], ...
    'LineWidth',0.8, 'FitBoxToText','on');

set(fig1,'PaperUnits','centimeters','PaperSize',[16 11]);

%% 4. SOLITON PARAMETERS ───────────────────────────────────────────────────
w = 1; m = 0.5; b = 1; a = 1; r = 0.5; c = -0.5;
d = 0.1; g = 0.2; e1 = 1; s1 = 1; l = 0; rho = 0.01;

n_core_min = min(unique_ncore);
n_eff_min  = clean_neff(unique_ncore == n_core_min);
n_core_max = max(unique_ncore);
n_eff_max  = clean_neff(unique_ncore == n_core_max);
n_core_mid = unique_ncore(round(end/2));
n_eff_mid  = clean_neff(round(end/2));

%% 5. SOLITON FUNCTION ────────────────────────────────────────────────────
function L = compute_soliton(x, t, n_core, n_eff, params)
    w = params.w; m = params.m; b = params.b; a = params.a;
    r = params.r; c = params.c; e1 = params.e1; s1 = params.s1;
    l = params.l; rho = params.rho; d = params.d; g = params.g;

    m_mod        = m * (1 - 0.5*(n_core - 1.33)/0.05);
    a_mod        = a * (1 + 2*(n_eff  - 1.33)/0.01);
    w_mod        = w * (1 - 1.5*(n_core - 1.33)/0.05);
    phase_factor = n_eff * 100;

    [X, T]        = meshgrid(x, t);
    theta         = w_mod * (X - m_mod * T);
    sh            = sinh(b*theta);
    ch            = cosh(b*theta);
    numerator     = 4 * a_mod * b * r * w_mod^2 * s1 * e1 * (sh - ch);
    inner_bracket = 1 + 2*ch.*sh - 2*ch.^2;
    denominator   = c*s1^2 + 4*a_mod^2*r^2*w_mod^2*e1^2*inner_bracket;
    phi           = m_mod*X + (T/2)*(b^2*w_mod^2 - m_mod^2 - 2*d*g) ...
                    + l + rho*params.B_t + phase_factor*X;
    L = (numerator ./ denominator) .* exp(1i*phi) + rho*randn(size(X));
end

%% 6. NUMERICAL GRID ──────────────────────────────────────────────────────
x  = linspace(-10, 10, 400);
t  = linspace(0, 10, 300);
[X, T] = meshgrid(x, t);
rng(42);
dt  = t(2) - t(1);
B_t = cumsum(randn(size(T)), 2) * sqrt(dt);
B_t = [zeros(size(B_t(:,1))), B_t(:,1:end-1)];

params.w=w; params.m=m; params.b=b; params.a=a; params.r=r; params.c=c;
params.e1=e1; params.s1=s1; params.d=d; params.g=g; params.l=l;
params.rho=rho; params.B_t=B_t;

%% 7. COMPUTE SOLITONS ────────────────────────────────────────────────────
fprintf('Computing solitons...\n');
L_low  = compute_soliton(x,t,n_core_min,n_eff_min,params); I_low  = abs(L_low).^2;
L_mid  = compute_soliton(x,t,n_core_mid,n_eff_mid,params); I_mid  = abs(L_mid).^2;
L_high = compute_soliton(x,t,n_core_max,n_eff_max,params); I_high = abs(L_high).^2;

%% ═══════════════════════════════════════════════════════════════════════
%  FIGURES 2–4 — WATERFALL PLOTS  (perceptually uniform colormap)
% ═══════════════════════════════════════════════════════════════════════
wf_titles = { ...
    sprintf('Low Glucose Concentration  (\\itn\\rm_{core} = %.6f)', n_core_min), ...
    sprintf('Moderate Glucose Concentration  (\\itn\\rm_{core} = %.6f)', n_core_mid), ...
    sprintf('High Glucose Concentration  (\\itn\\rm_{core} = %.6f)', n_core_max)};
wf_data   = {I_low, I_mid, I_high};
wf_pos    = [50 50; 100 50; 150 50];   % rows: [x y] per figure

% Build a smooth blue-orange-red colormap (better than jet for print)
cmap_wf = [linspace(0.05,0.95,256)', linspace(0.15,0.55,256)', linspace(0.55,0.05,256)'];

for k = 1:3
    fwf = figure('Units','centimeters','Position',[wf_pos(k,:) 18 13],'Color','w');
    Idata = wf_data{k};
    waterfall(X(1:25:end,:), T(1:25:end,:), Idata(1:25:end,:));
    colormap(cmap_wf);

    ax = gca;
    ax.Color       = [0.97 0.97 0.97];
    ax.GridAlpha   = 0.25;
    ax.GridColor   = [0.5 0.5 0.5];
    ax.XColor      = [0.2 0.2 0.2];
    ax.YColor      = [0.2 0.2 0.2];
    ax.ZColor      = [0.2 0.2 0.2];
    ax.TickDir     = 'out';
    ax.BoxStyle    = 'back';
    box on; grid on;
    view(40, 28);

    xlabel('Spatial Coordinate, \itx',    'FontSize',13);
    ylabel('Time, \itt',                  'FontSize',13);
    zlabel('Intensity |L_1|^2',           'FontSize',13);
    title(['Soliton Dynamics at ' wf_titles{k}], 'FontSize',13,'FontWeight','bold');

    cb = colorbar;
    cb.Label.String   = 'Intensity |L_1|^2';
    cb.Label.FontSize = 11;
    cb.TickDirection  = 'out';
    cb.LineWidth      = 0.8;

    set(fwf,'PaperUnits','centimeters','PaperSize',[18 13]);
end

%% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 5 — OSCILLOSCOPE 2D SUBPLOT  (3 stacked, black bg)
% ═══════════════════════════════════════════════════════════════════════
t_mid = round(length(t)/2);

col_low  = [0.18 0.80 0.44];   % phosphor green
col_mid  = [0.91 0.49 0.14];   % amber
col_high = [0.00 0.74 0.84];   % cyan
gc       = [0.28 0.28 0.28];   % grid color

sig_all  = {I_low(t_mid,:),  I_mid(t_mid,:),  I_high(t_mid,:)};
col_all  = {col_low,          col_mid,          col_high};
nc_all   = {n_core_min,       n_core_mid,       n_core_max};
lbl_all  = {'Low Glucose',    'Moderate Glucose','High Glucose'};

fig5 = figure('Units','centimeters','Position',[3 3 20 18],'Color','k');

for k = 1:3
    ax = subplot(3,1,k);
    sig = sig_all{k};
    pk  = max(sig);

    plot(x, sig, 'Color', col_all{k}, 'LineWidth', 2.0);
    hold on;

    % zero baseline
    plot(xlim,[0 0],'-','Color',[0.35 0.35 0.35],'LineWidth',0.6);

    set(ax, 'Color','k', ...
        'XColor',gc, 'YColor',gc, ...
        'FontSize',11, 'FontName','Times New Roman', ...
        'Box','on', 'LineWidth',0.8, ...
        'XGrid','on','YGrid','on', ...
        'GridAlpha',0.30,'GridLineStyle','-','GridColor',gc, ...
        'MinorGridAlpha',0.15,'XMinorGrid','on','YMinorGrid','on', ...
        'MinorGridLineStyle',':','MinorGridColor',gc, ...
        'TickDir','in','TickLength',[0.01 0.025]);

    xlim([min(x), max(x)]);
    ylim([0, pk*1.18]);

    % y-axis label
    ylabel('|L_1|^2', 'FontSize',12, 'Color','w');

    % Title bar inside plot (top-left)
    text(min(x)+0.4, pk*1.10, ...
        sprintf('%s    \\itn\\rm_{core} = %.6f', lbl_all{k}, nc_all{k}), ...
        'Color', col_all{k}, 'FontSize',11, 'FontWeight','bold', ...
        'FontName','Times New Roman');

    % Peak annotation arrow + value
    [~, pidx] = max(sig);
    annotation_x = x(pidx);
    annotation_y = pk;
    text(annotation_x + 0.3, pk*0.92, ...
        sprintf('\\uparrow  %.4f', pk), ...
        'Color',[1 1 1], 'FontSize',10, 'FontName','Times New Roman');

    % Propagation width (signal above 1% of peak)
    mask  = sig > pk * 0.01;
    if any(mask)
        xL = x(find(mask,1,'first'));
        xR = x(find(mask,1,'last'));
        pw = xR - xL;
        % draw horizontal span arrow
        mid_y = pk * 0.38;
        plot([xL xR],[mid_y mid_y],'-','Color',[0.55 0.55 0.55],'LineWidth',0.9);
        plot([xL xL],[mid_y-pk*0.04, mid_y+pk*0.04],'-','Color',[0.55 0.55 0.55],'LineWidth',0.9);
        plot([xR xR],[mid_y-pk*0.04, mid_y+pk*0.04],'-','Color',[0.55 0.55 0.55],'LineWidth',0.9);
        text((xL+xR)/2, mid_y + pk*0.07, sprintf('%.2f units', pw), ...
            'Color',[0.60 0.60 0.60],'FontSize',9,'HorizontalAlignment','center', ...
            'FontName','Times New Roman');
    end

    if k == 3
        xlabel('Spatial Coordinate, \itx', 'FontSize',12, 'Color','w');
    end
end

% Super-title
annotation('textbox',[0 0.965 1 0.035], ...
    'String','Oscilloscope View: Soliton Intensity  |L_1|^2  at Three Glucose Concentrations', ...
    'FontSize',12,'FontWeight','bold','Color','w', ...
    'HorizontalAlignment','center','EdgeColor','none','BackgroundColor','none', ...
    'FontName','Times New Roman');

set(fig5,'PaperUnits','centimeters','PaperSize',[20 18]);

%% ═══════════════════════════════════════════════════════════════════════
%  FIGURE 6 — SENSITIVITY vs SOLITON PROPERTIES  (dual panel)
% ═══════════════════════════════════════════════════════════════════════
peak_low   = max(I_low(:));
peak_mid   = max(I_mid(:));
peak_high  = max(I_high(:));
dx         = x(2) - x(1);
width_low  = sum(I_low(t_mid,:)  > peak_low/2)  * dx;
width_mid  = sum(I_mid(t_mid,:)  > peak_mid/2)  * dx;
width_high = sum(I_high(t_mid,:) > peak_high/2) * dx;

nc3    = [n_core_min, n_core_mid, n_core_max];
peaks3 = [peak_low,   peak_mid,   peak_high];
wids3  = [width_low,  width_mid,  width_high];

% Marker styles
ms = 9;
col_pts = [0.80 0.10 0.10;   % low  — red
           0.08 0.35 0.72;   % mid  — blue
           0.10 0.60 0.30];  % high — green

fig6 = figure('Units','centimeters','Position',[5 5 18 8],'Color','w');

% Panel (a)
ax6a = subplot(1,2,1);
hold on;
for k = 1:3
    plot(nc3(k), peaks3(k), 'o', ...
        'Color', col_pts(k,:), 'MarkerFaceColor', col_pts(k,:), ...
        'MarkerEdgeColor','w', 'MarkerSize', ms, 'LineWidth',1);
end
plot(nc3, peaks3, '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1.4);
for k = 1:3
    plot(nc3(k), peaks3(k), 'o', ...
        'Color', col_pts(k,:), 'MarkerFaceColor', col_pts(k,:), ...
        'MarkerEdgeColor','w', 'MarkerSize', ms, 'LineWidth',1);
end
grid on; box on;
ax6a.GridLineStyle  = ':';
ax6a.GridAlpha      = 0.40;
ax6a.TickDir        = 'out';
ax6a.XMinorTick     = 'on';
ax6a.YMinorTick     = 'on';
xlabel('Core Refractive Index, \itn\rm_{core}', 'FontSize',12);
ylabel('Peak Intensity  |L_1|^2_{max}',         'FontSize',12);
title('(a)  Peak Intensity vs Glucose Level',   'FontSize',12,'FontWeight','bold');
legend({'Low','Moderate','High'}, 'Location','best','FontSize',10,'Box','off');

% Panel (b)
ax6b = subplot(1,2,2);
hold on;
for k = 1:3
    plot(nc3(k), wids3(k), 's', ...
        'Color', col_pts(k,:), 'MarkerFaceColor', col_pts(k,:), ...
        'MarkerEdgeColor','w', 'MarkerSize', ms, 'LineWidth',1);
end
plot(nc3, wids3, '-', 'Color',[0.4 0.4 0.4], 'LineWidth',1.4);
for k = 1:3
    plot(nc3(k), wids3(k), 's', ...
        'Color', col_pts(k,:), 'MarkerFaceColor', col_pts(k,:), ...
        'MarkerEdgeColor','w', 'MarkerSize', ms, 'LineWidth',1);
end
grid on; box on;
ax6b.GridLineStyle  = ':';
ax6b.GridAlpha      = 0.40;
ax6b.TickDir        = 'out';
ax6b.XMinorTick     = 'on';
ax6b.YMinorTick     = 'on';
xlabel('Core Refractive Index, \itn\rm_{core}', 'FontSize',12);
ylabel('Pulse Width FWHM (spatial units)',       'FontSize',12);
title('(b)  Pulse Width vs Glucose Level',      'FontSize',12,'FontWeight','bold');
legend({'Low','Moderate','High'}, 'Location','best','FontSize',10,'Box','off');

sgtitle('Glucose Concentration Effects on Soliton Properties', ...
    'FontSize',13,'FontWeight','bold');

set(fig6,'PaperUnits','centimeters','PaperSize',[18 8]);


%% 11. AUTO-SAVE ALL FIGURES (600 dpi) ───────────────────────────────────
save_dir = fullfile(fileparts(mfilename('fullpath')), 'PCF_Figures');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end


% Grab figures in creation order
all_fig_h = flipud(findall(0,'Type','figure'));   % oldest first
fig_names = { ...
    'Fig1_Sensitivity_Analysis', ...
    'Fig2_Waterfall_Low_Glucose', ...
    'Fig3_Waterfall_Mid_Glucose', ...
    'Fig4_Waterfall_High_Glucose', ...
    'Fig5_Oscilloscope_2D', ...
    'Fig6_Soliton_Properties'};

% Figs 1,5,6 are 2D — vector PDF sharp; Figs 2-4 are 3D waterfall — image PDF
vector_figs = [1, 5, 6];

fprintf('\nSaving figures — 600 dpi PNG / TIFF / PDF...\n');
for k = 1 : min(numel(all_fig_h), numel(fig_names))
    fh   = all_fig_h(k);
    fout = fullfile(save_dir, fig_names{k});

    % PNG 600 dpi
    exportgraphics(fh, [fout '.png'], ...
        'Resolution', 600, 'BackgroundColor', 'current');

    % TIFF 600 dpi (Elsevier preferred)
    exportgraphics(fh, [fout '.tif'], ...
        'Resolution', 600, 'BackgroundColor', 'current');

    % PDF: vector for 2D, high-res image for 3D (no warning)
    if ismember(k, vector_figs)
        exportgraphics(fh, [fout '.pdf'], ...
            'ContentType', 'vector', 'BackgroundColor', 'current');
    else
        exportgraphics(fh, [fout '.pdf'], ...
            'ContentType', 'image', 'Resolution', 600, ...
            'BackgroundColor', 'current');
    end

    fprintf('  [%d/%d] %s\n', k, numel(fig_names), fig_names{k});
end
fprintf('Done. Figures saved to:\n  %s\n', save_dir);

%% 12. OUTPUT SUMMARY ─────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('Q1 JOURNAL PAPER - RESULTS SUMMARY\n');
fprintf('========================================\n');
fprintf('Sensitivity (S) = %.4f RIU^{-1}\n', Sensitivity);
fprintf('R-squared (R^2) = %.6f\n', R_sq);
fprintf('\nGlucose Effect on Soliton:\n');
fprintf('Low  Glucose (n_core=%.5f): Peak=%.4f, FWHM=%.4f\n', n_core_min, peak_low,  width_low);
fprintf('Mid  Glucose (n_core=%.5f): Peak=%.4f, FWHM=%.4f\n', n_core_mid, peak_mid,  width_mid);
fprintf('High Glucose (n_core=%.5f): Peak=%.4f, FWHM=%.4f\n', n_core_max, peak_high, width_high);
fprintf('========================================\n');
