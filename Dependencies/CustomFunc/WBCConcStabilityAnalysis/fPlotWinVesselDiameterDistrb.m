function [Plot_WinSz_Out] = fPlotWinVesselDiameterDistrb(GroupID_AllVidData_WBCConcStability_Tbl_In,UI_In)

Diameter_In = GroupID_AllVidData_WBCConcStability_Tbl_In.WinVesselDiameter_um;

% = Histogram [subplot2 REMOVED]
FigHandle = figure('Position',[30 500 900 300],'visible','on');
tlo = tiledlayout(1,1,'TileSpacing','Compact','Padding','tight'); % tile layout; replace subplot()

ax = nexttile(tlo); 
hold(ax, 'on');
grid on;
histogram(Diameter_In,[0.5*floor(min(Diameter_In,[],'all')/0.5):0.5:0.5*ceil(max(Diameter_In,[],'all')/0.5)],'FaceColor',[0.9290 0.6940 0.1250]);
xlim([0.5*floor(min(Diameter_In,[],'all')/0.5) 0.5*ceil(max(Diameter_In,[],'all')/0.5)]);
xticks([0.5*floor(min(Diameter_In,[],'all')/0.5):0.5:0.5*ceil(max(Diameter_In,[],'all')/0.5)]);
xlabel('WinVesselDiameter (um)','interpreter','none');
ylabel('HistCt','interpreter','none');
title({strcat('# CellCt Window Vessel Diameter values:',32,num2str(numel(Diameter_In))),...
    strcat('CellCt Window Vessel Diameter HistEdge, User-Input (um):',32,strjoin(arrayfun(@num2str,UI_In.VesselDiameterHistEdge_um,'UniformOutput',false)))},...
    'fontsize',6,'fontweight','normal','interpreter','none');
fontsize(ax,12,"points");

% = Export
Plot_WinSz_Out = export_fig('-r150');

close all;
clearvars tlo ax ax2 cdf_handle FigHandle;
