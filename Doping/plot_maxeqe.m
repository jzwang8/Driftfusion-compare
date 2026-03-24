% Data ---------------------------------------------------------------
eqe = [	84.90719328	115.3121491; 84.93291074	115.4067051; 84.94822135	115.4628533];  % rows: Qf, cols: [eta=0, eta=2]
% % voc
% eqe = 
% % pce
% 
% % jsc


fixedCharge = {'-1×10^{19}', '0', '+1×10^{19}'};
x = 1:numel(fixedCharge);

c0 = [0.33 0.51 0.66];   % blue   (η_SF=0)
c2 = [1.00 0.55 0.07];   % orange (η_SF=2)

% Figure --------------------------------------------------------------
figure('Color','w','Position',[120 120 560 320]); hold on;

% connectors
for i = 1:numel(x)
    plot([x(i) x(i)], [eqe(i,1) eqe(i,2)], '-', 'Color',[0.75 0.75 0.75], 'LineWidth',2);
end

% endpoints
s0 = scatter(x, eqe(:,1), 64, c0, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
s2 = scatter(x, eqe(:,2), 64, c2, 'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);

% axes & labels
ylim([70 140]); xlim([0.5 3.5]);
set(gca,'XTick',x,'XTickLabel',fixedCharge,'Box','on','LineWidth',1.2,'FontSize',14);
grid off
xlabel('Fixed charge density (cm^{-3})','FontSize',14);
ylabel('EQE (%)','FontSize',14);

% legend (below)
legend([s0 s2], {'\eta_{SF} = 0','\eta_{SF} = 2'}, ...
       'Location','northwest','Orientation','vertical', ...
       'FontSize',12,'Box','on','Interpreter','tex');

% value labels: η=0 below, η=2 above (with dynamic offset)
yl = ylim; dy = 0.02*(yl(2)-yl(1));   % 2% of y-range
for i = 1:numel(x)
    % ensure spacing even if η2 < η0 (robust)
    yLow  = min(eqe(i,1), eqe(i,2));
    yHigh = max(eqe(i,1), eqe(i,2));
    % place η=0 below
    text(x(i), eqe(i,1) - dy, sprintf('%.1f', eqe(i,1)), ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',12);
    % place η=2 above
    text(x(i), eqe(i,2) + dy, sprintf('%.1f', eqe(i,2)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12);
end

% % -------- Max-EQE data ---------------------------------------------
% % rows → fixed-charge density   (-1e19, 0, +1e19  cm-3)
% % cols → SF efficiency          (0 and 2)
% eqe = [ ...
%     76.41   107.927 ;   % −1e19
%     78.5137   118.169 ;   % 0
%     78.734   119.242 ];  % +1e19
% 
% fixedCharge = {'-1×10^{19}', '0', '+1×10^{19}'};
% sfSeries    = {'\eta_{SF} = 0', '\eta_{SF} = 2'};          % legend follows matrix columns
% 
% % -------- Clustered bar chart --------------------------------------
% figure('Color','w')
% b   = bar(eqe,'grouped');
% 
% ax  = gca;
% ax.FontSize = 18;                           % bigger tick‐label font
% ax.XTickLabel = fixedCharge;
% ax.XTickLabelRotation = 0;                  % keep tick labels horizontal
% ylim([0,130])
% ylabel('Max EQE (%)','FontSize',18)
% xlabel('Fixed charge density (cm^{-3})','FontSize',18)
% 
% legend(sfSeries,'Location','northwest','FontSize',16)
% 
% % -------- Annotate each bar with its value (horizontal) ------------
% for i = 1:numel(b)                          % loop over the two SF groups
%     x = b(i).XEndPoints;
%     y = b(i).YEndPoints;
%     labels = compose('%.3f', b(i).YData);   % keeps three decimals
%     text(x, y + 1, labels, ...              % +1 sets a small gap above bar
%          'HorizontalAlignment','center', ...
%          'VerticalAlignment','bottom', ...
%          'FontSize',18, ...
%          'Rotation',0);                     % horizontal text
% end