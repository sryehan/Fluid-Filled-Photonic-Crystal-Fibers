%% ============================================================
%% Phase Portrait Analysis
%% w^2*M'' - (m^2 + 2(n+dg))*M - 2c*M^3 = 0
%% System: M' = V,  V' = (1/w^2)*(alpha*M + beta*M^3)
%% Hamiltonian: H = 0.5*V^2 - (1/w^2)*(alpha/2*M^2 + beta/4*M^4)
%% ============================================================
clear; clc; close all;

%% ---- Scenario Definitions ----
% Each row: [w, m, n, d, g, c]
S(1).w=1.0; S(1).m=0.5; S(1).n=-0.8; S(1).d=0.1; S(1).g=0.2; S(1).c=-0.5;
S(1).name='S1: Bright Soliton (a<0,c<0) [BEST]';
S(1).col=[0 0.45 0.74];

S(2).w=1.0; S(2).m=1.0; S(2).n=0.3;  S(2).d=0.2; S(2).g=0.5; S(2).c=-1.0;
S(2).name='S2: Dark/Kink Soliton (a>0,c<0) [GOOD]';
S(2).col=[0.17 0.63 0.17];

S(3).w=1.0; S(3).m=1.0; S(3).n=-0.5; S(3).d=0.0; S(3).g=0.0; S(3).c=-0.5;
S(3).name='S3: Marginal Pure Cubic (a=0,c<0)';
S(3).col=[0.93 0.69 0.13];

S(4).w=1.0; S(4).m=1.5; S(4).n=0.5;  S(4).d=0.3; S(4).g=0.4; S(4).c=0.5;
S(4).name='S4: Explosive (a>0,c>0) [WORST]';
S(4).col=[0.85 0.11 0.11];

S(5).w=1.0; S(5).m=0.3; S(5).n=-1.5; S(5).d=0.1; S(5).g=0.1; S(5).c=2.0;
S(5).name='S5: Defocusing (a<0,c>0) [BAD]';
S(5).col=[0.64 0.08 0.18];

S(6).w=2.5; S(6).m=0.5; S(6).n=-0.9; S(6).d=0.1; S(6).g=0.2; S(6).c=-0.3;
S(6).name='S6: Wide Soliton Large w (a<0,c<0) [BEST+]';
S(6).col=[0.49 0.18 0.56];

%% ============================================================
%% LOOP: One figure per scenario
%% ============================================================
for k = 1:6
    sc   = S(k);
    alph = sc.m^2 + 2*(sc.n + sc.d*sc.g);   % alpha
    beta = 2*sc.c;                            % beta
    cf   = 1/sc.w^2;                          % 1/w^2

    % --- Critical points ---
    cps = 0;
    if beta ~= 0 && (-alph/beta) > 0
        Mex = sqrt(-alph/beta);
        cps = [0, Mex, -Mex];
    end

    % --- Axis range ---
    Mrng = max(abs(cps))*2.8;
    if Mrng < 2.5; Mrng = 2.5; end
    Vrng = Mrng*1.1;

    % --- Meshgrid for quiver ---
    [MM,VV] = meshgrid(linspace(-Mrng,Mrng,26), ...
                       linspace(-Vrng,Vrng,26));
    dM =  VV;
    dV =  cf*(alph*MM + beta*MM.^3);
    mag = sqrt(dM.^2+dV.^2); mag(mag==0)=1;
    dMn = dM./mag;  dVn = dV./mag;

    % --- Hamiltonian for contours ---
    Mf = linspace(-Mrng*1.05, Mrng*1.05, 600);
    Vf = linspace(-Vrng*1.05, Vrng*1.05, 600);
    [Mfg,Vfg] = meshgrid(Mf,Vf);
    H = 0.5*Vfg.^2 - cf*(alph/2*Mfg.^2 + beta/4*Mfg.^4);

    Hmin = prctile(H(:),2);
    Hmax = prctile(H(:),98);
    lvls = linspace(Hmin, Hmax, 45);

    % Separatrix levels
    sep_lvls = [];
    for cp = cps
        if abs(cp) > 1e-9
            Hcp = -cf*(alph/2*cp^2 + beta/4*cp^4);
            if Hcp > Hmin && Hcp < Hmax
                sep_lvls(end+1) = Hcp;
            end
        end
    end

    % --- FIGURE ---
    figure('Position',[80 80 880 720],'Color','w');

    % Trajectory contours
    contour(Mfg, Vfg, H, lvls, 'LineWidth', 1.1);
    colormap(jet(256));
    caxis([Hmin Hmax]);
    hold on;

    % Separatrix (dashed black)
    for sl = sep_lvls
        contour(Mfg, Vfg, H, [sl sl], 'k--', 'LineWidth', 2.5);
    end

    % Quiver
    quiver(MM, VV, dMn, dVn, 0.5, 'Color', [0.45 0.45 0.45], ...
           'LineWidth', 0.7, 'MaxHeadSize', 0.5);

    % Critical points + labels
    for cp = cps
        J22 = cf*(alph + 3*beta*cp^2);
        if J22 > 1e-9
            cptype = 'SADDLE'; mk='rs';
        elseif J22 < -1e-9
            cptype = 'CENTER'; mk='bo';
        else
            cptype = 'DEGEN.'; mk='k^';
        end
        plot(cp, 0, mk, 'MarkerSize',11, ...
             'MarkerFaceColor','yellow','LineWidth',2);
        text(cp+0.05*Mrng, 0.1*Vrng, ...
             sprintf('(%.2f,0)\n%s',cp,cptype), ...
             'FontSize',8,'FontWeight','bold');
    end

    % Axis lines
    plot([-Mrng Mrng],[0 0],'k-','LineWidth',0.9);
    plot([0 0],[-Vrng Vrng],'k-','LineWidth',0.9);

    xlabel('\mathcal{M}  (amplitude)','FontSize',13,'FontWeight','bold');
    ylabel('\mathcal{M}''  =  V  (slope)','FontSize',13,'FontWeight','bold');
    xlim([-Mrng Mrng]); ylim([-Vrng Vrng]);
    grid on; box on;
    colorbar('Location','EastOutside');

    % Title
    title({sc.name, ...
           sprintf('\\alpha=%.3f   \\beta=2c=%.3f   w=%.2f', alph,beta,sc.w)}, ...
          'FontSize',11,'FontWeight','bold','Color',sc.col);

    % Bottom annotation
    verdict = get_verdict(k);
    dim = [0.10 0.01 0.82 0.07];
    annotation('textbox',dim,'String', ...
        {sprintf('Params: m=%.2f  n=%.2f  d=%.2f  g=%.2f  c=%.2f  w=%.2f', ...
                  sc.m,sc.n,sc.d,sc.g,sc.c,sc.w), verdict}, ...
        'FontSize',9,'EdgeColor',sc.col,'BackgroundColor',[0.97 0.97 1], ...
        'LineWidth',1.8,'FitBoxToText','off');

    hold off;
    drawnow;
end

%% ============================================================
%% PARAMETER REGIME MAP  (alpha-c plane)
%% ============================================================
figure('Position',[80 80 950 760],'Color','w');

av = linspace(-3,3,500);
cv = linspace(-3,3,500);
[AA,CC] = meshgrid(av,cv);

region = zeros(size(AA));
region((AA<0)&(CC<0)) = 1;   % bright soliton
region((AA>0)&(CC<0)) = 2;   % dark soliton
region((AA<0)&(CC>0)) = 3;   % defocusing
region((AA>0)&(CC>0)) = 4;   % explosive

cmap4 = [0.20 0.60 1.00;   % blue   - bright
         0.20 0.80 0.20;   % green  - dark
         1.00 0.80 0.00;   % yellow - defocusing
         0.90 0.20 0.20];  % red    - explosive
colormap(cmap4);

imagesc(av, cv, region);
set(gca,'YDir','normal');
caxis([0.5 4.5]);
cb = colorbar('Ticks',[1 2 3 4], ...
    'TickLabels',{'1: Bright Soliton  [BEST]', ...
                  '2: Dark/Kink Soliton [GOOD]', ...
                  '3: Defocusing  [POOR]', ...
                  '4: Explosive  [WORST]'});
cb.FontSize = 10;

hold on;
plot([-3 3],[0 0],'k-','LineWidth',2);
plot([0 0],[-3 3],'k-','LineWidth',2);

% Mark scenarios on map
for k = 1:6
    sc   = S(k);
    ak   = sc.m^2 + 2*(sc.n + sc.d*sc.g);
    plot(ak, sc.c, 'w*','MarkerSize',16,'LineWidth',2);
    text(ak+0.12, sc.c+0.12, sprintf('S%d',k), ...
         'Color','white','FontWeight','bold','FontSize',11, ...
         'BackgroundColor','black');
end

% Region labels
text(-1.7,-2.1,{'BRIGHT SOLITON','✓ USE THIS','Closed orbits'}, ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold', ...
    'Color','white','BackgroundColor',[0 0.3 0.8]);
text(1.6,-2.1,{'DARK/KINK','✓ VALID','Heteroclinic orbit'}, ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold', ...
    'Color','white','BackgroundColor',[0.1 0.5 0.1]);
text(-1.7,2.1,{'DEFOCUSING','✗ AVOID','Small amp only'}, ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold', ...
    'Color','black','BackgroundColor',[0.9 0.75 0]);
text(1.6,2.1,{'EXPLOSIVE','✗ WORST','No soliton!'}, ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold', ...
    'Color','white','BackgroundColor',[0.75 0.1 0.1]);

xlabel('\alpha = m^2 + 2(n+dg)   [Linear coefficient]', ...
       'FontSize',13,'FontWeight','bold');
ylabel('c   [Nonlinear coefficient]','FontSize',13,'FontWeight','bold');
title({'Parameter Space: Physical Regions', ...
       'w^2 M'''' - \alpha M - 2c M^3 = 0', ...
       'Choose parameters from BLUE region (S1,S2,S6) for soliton solutions'}, ...
      'FontSize',12,'FontWeight','bold');

xlim([-3 3]); ylim([-3 3]);
grid on;
hold off;
drawnow;

%% ============================================================
%% CONSOLE: Parameter Guide
%% ============================================================
fprintf('\n=========================================================\n');
fprintf('   PARAMETER GUIDE\n');
fprintf('   alpha = m^2+2(n+dg),   beta = 2c\n');
fprintf('=========================================================\n');
fprintf('%-20s %-10s %-10s %s\n','REGIME','alpha','c','Physical result');
fprintf('%s\n',repmat('-',1,60));
fprintf('%-20s %-10s %-10s %s\n','Bright soliton','< 0','< 0','BEST: localized solution');
fprintf('%-20s %-10s %-10s %s\n','Dark soliton','> 0','< 0','GOOD: kink/dark soliton');
fprintf('%-20s %-10s %-10s %s\n','Defocusing','< 0','> 0','POOR: small amplitude only');
fprintf('%-20s %-10s %-10s %s\n','Explosive','> 0','> 0','WORST: no soliton');
fprintf('\nSoliton Amplitude : A = sqrt(|alpha|/(2|c|))\n');
fprintf('Soliton Width     : proportional to  w/sqrt(|alpha|)\n');
fprintf('\nRecommended (bright soliton):\n');
fprintf('  alpha in [-3.0, -0.1]  =>  n + d*g < -m^2/2\n');
fprintf('  c     in [-2.0, -0.1]\n');
fprintf('  w     in [0.5,  3.0]\n');
fprintf('=========================================================\n\n');

%% ============================================================
%% Helper function
%% ============================================================
function v = get_verdict(k)
    verdicts = {
        'BEST: alpha<0, c<0 => bounded closed orbits => bright soliton solution', ...
        'GOOD: alpha>0, c<0 => saddle + heteroclinic orbit => dark/kink soliton', ...
        'MARGINAL: alpha~0, c<0 => nonlinear center, periodic orbits only', ...
        'WORST A: alpha>0, c>0 => all orbits diverge, no physical soliton', ...
        'BAD: alpha<0, c>0 => defocusing, solution escapes for large M', ...
        'BEST+: large w, alpha<0, c<0 => wide smooth soliton, stable orbits'};
    v = verdicts{k};
end