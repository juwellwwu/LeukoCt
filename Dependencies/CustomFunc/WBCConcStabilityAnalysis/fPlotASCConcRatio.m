function [Plot_ConcRatio_Cell_Out] = fPlotASCConcRatio(GroupID_AllVidData_WBCConcStability_Tbl_In,UI_In)

fprintf('Plotting boxplot of WBC concentration ratio (AUTO / REF), grouped by cell count window diameter ...\n');

Plot_ConcRatio_Cell_Out = cell(1,1,2);


%%

VesselDiameterHistEdge_In = UI_In.VesselDiameterHistEdge_um;


%% Plot 01: Boxplot

% === Data: Plot 1 (for Classic MATLAB boxplot), Plot 3 (gscatter,
% WBCRatio), Plot 4 (gscatter, MaxWBCConcRatioDiff)
bpData_Tbl = GroupID_AllVidData_WBCConcStability_Tbl_In;
bpData_Tbl(bpData_Tbl.ASCImgTimeStart_s<0,:) = []; 
bpData_Tbl(bpData_Tbl.WinVesselDiameter_um<min(VesselDiameterHistEdge_In,[],'all'),:) = [];
bpData_Tbl(bpData_Tbl.WinVesselDiameter_um>max(VesselDiameterHistEdge_In,[],'all'),:) = [];

if UI_In.ASCOption==1
    bpData = bpData_Tbl.FinalIncrVolBinWBCConc_KuL./bpData_Tbl.REFWBCConc_KuL;
elseif UI_In.ASCOption==2
    bpData = bpData_Tbl.FinalCumVolBinWBCConc_KuL./bpData_Tbl.REFWBCConc_KuL;
end

% = Divide boxplot by Window Vessel Diameter (um) 
VesselDiameterHistBinCtr = 0.5*movsum(VesselDiameterHistEdge_In,2,'Endpoints','discard');
[VesselDiameterHistCt,~,VesselDiameterHistBin] = histcounts(bpData_Tbl.WinVesselDiameter_um,VesselDiameterHistEdge_In);
[grp_VesselDiameterHistBin, grpID_VesselDiameterHistBin] = findgroups(VesselDiameterHistBin);
cmap = colormap(jet(numel(VesselDiameterHistEdge_In)-1)); % Rainbow; Small diameter = Blue; Large diameter = Red

% = Statistics by vessel diameter group (to be placed in title)
Ratio_Mean = splitapply(@mean,bpData,grp_VesselDiameterHistBin);
Ratio_Mean = horzcat(VesselDiameterHistBinCtr(grpID_VesselDiameterHistBin),Ratio_Mean);

Ratio_Stdev = splitapply(@std,bpData,grp_VesselDiameterHistBin);
Ratio_Stdev = horzcat(VesselDiameterHistBinCtr(grpID_VesselDiameterHistBin),Ratio_Stdev);

Ratio_Median = splitapply(@median,bpData,grp_VesselDiameterHistBin);
Ratio_Median = horzcat(VesselDiameterHistBinCtr(grpID_VesselDiameterHistBin),Ratio_Median);

Ratio_25prctile = splitapply(@(x) prctile(x,25),bpData,grp_VesselDiameterHistBin);
Ratio_25prctile = horzcat(VesselDiameterHistBinCtr(grpID_VesselDiameterHistBin),Ratio_25prctile);

Ratio_75prctile = splitapply(@(x) prctile(x,75),bpData,grp_VesselDiameterHistBin);
Ratio_75prctile = horzcat(VesselDiameterHistBinCtr(grpID_VesselDiameterHistBin),Ratio_75prctile);

% === Data: Plot 2 (for daboxplot)
% function daboxplot() from FileExchange:
% https://www.mathworks.com/matlabcentral/fileexchange/136524-daviolinplot-beautiful-violin-and-raincloud-plots

bpData2_Cell = sortrows(horzcat(bpData,grp_VesselDiameterHistBin),2,'ascend'); 
bpData2_Cell(:,2) = []; 
bpData2_Cell = permute(mat2cell(bpData2_Cell,histcounts(grp_VesselDiameterHistBin)),[2 1]);

% = Plot
FigHandle = figure('Position',[20 250 450 550],'visible','on');
tlo = tiledlayout(1,1,'TileSpacing','Compact','Padding','tight'); 

% % % ax = nexttile(tlo); 
% % % hold(ax, 'on');
% % % grid on;
% % % 
% % % boxplot(ax, bpData,grp_VesselDiameterHistBin,'labels',{num2str(VesselDiameterHistBinCtr(grpID_VesselDiameterHistBin))});
% % % xlabel('WinVesselDiameter_HistBinCtr (um)','interpreter','none')
% % % if UI_In.ASCOption==1
% % %     ylabel('IncrWBCConc_FinalVolBin/REFConc','interpreter','none');
% % % elseif UI_In.ASCOption==2
% % %     ylabel('CumWBCConc_FinalVolBin/REFConc','interpreter','none');
% % % end
% % % ylim([0.5*floor(min(bpData,[],1)/0.5) 0.5*ceil(max(bpData,[],1)/0.5)]);
% % % title({strcat('CellCt Window Vessel Diameter HistEdge (um):',32,strjoin(arrayfun(@num2str,VesselDiameterHistEdge_In,'UniformOutput',false))),...
% % %     strcat('CellCt Window Vessel Diameter HistBinCtr (um):',32,strjoin(arrayfun(@num2str,VesselDiameterHistBinCtr,'UniformOutput',false))),...
% % %     strcat('# ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,VesselDiameterHistCt,'UniformOutput',false))),...
% % %     strcat('Ratio Mean, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_Mean(:,2),'UniformOutput',false))),...
% % %     strcat('Ratio Stdev, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_Stdev(:,2),'UniformOutput',false))),...
% % %     strcat('Ratio Median, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_Median(:,2),'UniformOutput',false))),...
% % %     strcat('Ratio 25th prctile, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_25prctile(:,2),'UniformOutput',false))),...
% % %     strcat('Ratio 75th prctile, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_75prctile(:,2),'UniformOutput',false))),...
% % %     },...
% % %     'fontsize',6,'fontweight','normal','interpreter','none');

% % % ax2 = nexttile(tlo); 
% % % hold(ax2, 'on');
% % % grid on;
% % % 
% % % daboxplot(bpData2_Cell,...
% % %     'xtlabels',string((VesselDiameterHistBinCtr(grpID_VesselDiameterHistBin))'),... % conditions
% % %     'flipcolors',1,'colors',cmap,...
% % %     'fill',0,'boxwidth',1.2,'boxspacing',1.2,...
% % %     'scatter',2,'scattersize',20,'scatteralpha',0.8); % 'scatter',1 doesn't work
% % % 
% % % if UI_In.ASCOption==1
% % %     ylabel('IncrWBCConc_FinalVolBin/REFConc','interpreter','none');
% % % elseif UI_In.ASCOption==2
% % %     ylabel('CumWBCConc_FinalVolBin/REFConc','interpreter','none');
% % % end
% % % ylim([0.5*floor(min(bpData,[],1)/0.5) 0.5*ceil(max(bpData,[],1)/0.5)]);
% % % 
% % % title({strcat('CellCt Window Vessel Diameter HistEdge (um):',32,strjoin(arrayfun(@num2str,VesselDiameterHistEdge_In,'UniformOutput',false))),...
% % %     strcat('CellCt Window Vessel Diameter HistBinCtr (um):',32,strjoin(arrayfun(@num2str,VesselDiameterHistBinCtr,'UniformOutput',false))),...
% % %     strcat('# ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,VesselDiameterHistCt,'UniformOutput',false))),...
% % %     strcat('Ratio Mean, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_Mean(:,2),'UniformOutput',false))),...
% % %     strcat('Ratio Stdev, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_Stdev(:,2),'UniformOutput',false))),...
% % %     strcat('Ratio Median, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_Median(:,2),'UniformOutput',false))),...
% % %     strcat('Ratio 25th prctile, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_25prctile(:,2),'UniformOutput',false))),...
% % %     strcat('Ratio 75th prctile, ASC-Stable GroupIDs in HistBins:',32,strjoin(arrayfun(@num2str,Ratio_75prctile(:,2),'UniformOutput',false))),...
% % %     },...
% % %     'fontsize',6,'fontweight','normal','interpreter','none');

ax3 = nexttile(tlo); 
hold(ax3, 'on');
grid on;

gscatter(bpData_Tbl.WinVesselDiameter_um,bpData,grp_VesselDiameterHistBin,cmap);
labelpoints(bpData_Tbl.WinVesselDiameter_um,bpData,bpData_Tbl.GroupID,'FontSize',6);
xlabel('WinVesselDiameter (um)','interpreter','none')
if UI_In.ASCOption==1
    ylabel('IncrWBCConc_FinalVolBin/REFConc','interpreter','none');
elseif UI_In.ASCOption==2
    ylabel('CumWBCConc_FinalVolBin/REFConc','interpreter','none');
end
ylim([0.5*floor(min(bpData,[],1)/0.5) 0.5*ceil(max(bpData,[],1)/0.5)]);

% = Plot Output
Plot_ConcRatio_Cell_Out{:,:,1} = export_fig('-r150');

clearvars tlo ax*;


%% Plot 02 [Removed]

Plot_ConcRatio_Cell_Out{:,:,2} = [];

clearvars tlo ax* color_*;

%%

close all;
clearvars bpData* FigHandle ax* tlo VesselDiameterHistBinCtr VesselDiameterHistCt VesselDiameterHistBin grp_VesselDiameterHistBin grpID_VesselDiameterHistBin;
