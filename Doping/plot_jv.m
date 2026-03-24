%% Dumbbell plots: SF = 0 vs SF = 2, side-by-side Gaussian vs Exponential
close all; clearvars;

% ---- Data (order: Qf = [+1e19, -1e19, 0]) ----
% Voc (V)
Voc_G_SF0 = [0.6533,  0.64283, 0.65122];
Voc_G_SF2 = [0.65486, 0.64344, 0.65250];
Voc_E_SF0 = [0.65793, 0.65781, 0.65788];
Voc_E_SF2 = [0.65965, 0.65953, 0.65961];

% PCE (%)
PCE_G_SF0 = [20.21,   18.741,  19.933];
PCE_G_SF2 = [21.46,   19.2073, 21.0518];
PCE_E_SF0 = [20.7432, 20.7358, 20.7404];
PCE_E_SF2 = [22.2805, 22.2671, 22.2754];

% Jsc (mA cm^-2)
Jsc_G_SF0 = [36.9638, 34.9374, 36.6049];
Jsc_G_SF2 = [39.1711, 35.7695, 38.5679];
Jsc_E_SF0 = [37.6481, 37.6436, 37.6470];
Jsc_E_SF2 = [40.3135, 40.2993, 40.3086];

% ---- Plot settings ----
Qf_labels = {'+1e19','-1e19','0'};
x = 1:numel(Qf_labels);
dx = 0.2;                          % small offset to separate G vs E
xG = x - dx;                       % Gaussian x-positions
xE = x + dx;                       % Exponential x-positions
cG = [0.00 0.45 0.74];             % Gaussian color
cE = [0.85 0.33 0.10];             % Exponential color
ms = 50; lw = 2;

figure('Color','w','Units','inches','Position',[0.5 0.5 11 3.6]);
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- Helper to draw one panel ----
drawPanel = @(ax, yG0,yG2,yE0,yE2, ylab) ...
    localDumbbell(ax, xG,xE, yG0,yG2,yE0,yE2, Qf_labels, cG,cE, ms,lw, ylab);

% Voc
ax1 = nexttile; 
drawPanel(ax1, Voc_G_SF0, Voc_G_SF2, Voc_E_SF0, Voc_E_SF2, 'V_{OC} (V)');

% PCE
ax2 = nexttile; 
drawPanel(ax2, PCE_G_SF0, PCE_G_SF2, PCE_E_SF0, PCE_E_SF2, 'PCE (%)');

% Jsc
ax3 = nexttile; 
drawPanel(ax3, Jsc_G_SF0, Jsc_G_SF2, Jsc_E_SF0, Jsc_E_SF2, 'J_{SC} (mA cm^{-2})');

% title(t, 'Effect of singlet fission (SF = 0 \rightarrow 2) at three Q_f values; Gaussian vs Exponential','FontWeight','normal');

% ---- Legend (build once from dummy lines) ----
axes(ax1); hold on;
hG0 = plot(nan,nan,'-','Color',cG,'LineWidth',lw); % line style for G
hG2 = scatter(nan,nan,ms,cG,'filled','MarkerEdgeColor','k'); % filled endpoint
hE0 = plot(nan,nan,'-','Color',cE,'LineWidth',lw);
hE2 = scatter(nan,nan,ms,cE,'filled','MarkerEdgeColor','k');
legend([hG0 hE0], ...
       {'Gaussian (SF=0 \leftrightarrow SF=2)',...
        'Exponential (SF=0 \leftrightarrow SF=2)'},...
       'Location','southoutside','NumColumns',2,'Box','off');

% ===== Local function =====
function localDumbbell(ax, xG,xE, yG0,yG2,yE0,yE2, xlabels, cG,cE, ms,lw, ylab)
    axes(ax); cla; hold on;
    % Gaussian dumbbells
    for i = 1:numel(xG)
        plot([xG(i) xG(i)],[yG0(i) yG2(i)],'-','Color',cG,'LineWidth',lw);
        scatter(xG(i), yG0(i), ms, 'o', 'MarkerEdgeColor', cG, 'MarkerFaceColor','w');
        scatter(xG(i), yG2(i), ms, cG, 'filled', 'MarkerEdgeColor','k');
    end
    % Exponential dumbbells
    for i = 1:numel(xE)
        plot([xE(i) xE(i)],[yE0(i) yE2(i)],'-','Color',cE,'LineWidth',lw);
        scatter(xE(i), yE0(i), ms, 'o', 'MarkerEdgeColor', cE, 'MarkerFaceColor','w');
        scatter(xE(i), yE2(i), ms, cE, 'filled', 'MarkerEdgeColor','k');
    end
    % Cosmetics
    xlim([min([xG xE])-0.4, max([xG xE])+0.4]);
    set(ax,'XTick',1:numel(xlabels),'XTickLabel',xlabels,'FontSize',11,'LineWidth',1);
    xlabel('Q_f (cm^{-3})'); ylabel(ylab);
    grid on; box on;
end