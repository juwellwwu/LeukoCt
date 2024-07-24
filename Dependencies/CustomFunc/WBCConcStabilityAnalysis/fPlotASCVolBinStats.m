function [Plot_ASCVolBinStats_Out] = fPlotASCVolBinStats(GroupID_AllVidData_WBCConcStability_Tbl_In,UI_In)

fprintf('Plotting Statistics (boxplot, histogram) of Starting ImgTime (s) and CumFlowVol (nL) at which ASC-determined Stable WBC Conc begins ...\n');

%%

VesselDiameterHistEdge_In = UI_In.VesselDiameterHistEdge_um;

if UI_In.ASCOption == 1 
    ASCRatioDiffLim = UI_In.ASC.IncrVolBin.ASCRatioDiffLim; 
    ASC_ReportIncVol_Step_In = UI_In.ASC.IncrVolBin.ReportIncVol_Step;    
elseif UI_In.ASCOption == 2 
    ASCRatioDiffLim = UI_In.ASC.CumVolBin.ASCRatioDiffLim;
    ASC_ReportIncVol_Step_In = UI_In.ASC.CumVolBin.ReportIncVol_Step;
end


%%

% = Boxplot [subplots 3-5 REMOVED]
FigHandle = figure('Position',[100 250 640 750],'visible','on');
tlo = tiledlayout(1,2,'TileSpacing','Compact','Padding','tight'); % tile layout; replace subplot()

ax = nexttile(tlo); 
hold(ax, 'on');
grid on;
bpData = GroupID_AllVidData_WBCConcStability_Tbl_In.ASCImgTimeStart_s(GroupID_AllVidData_WBCConcStability_Tbl_In.ASCImgTimeStart_s>eps);
boxplot(ax, bpData);
xlabel('All_GroupID','interpreter','none')
xlim([0.5 1.5]);
if UI_In.ASCOption==1
    ylabel('ImgTimeStart_s (IncrWBCConc)','interpreter','none');
elseif UI_In.ASCOption==2
    ylabel('ImgTimeStart_s (CumWBCConc)','interpreter','none');
end
ylim([0 max(GroupID_AllVidData_WBCConcStability_Tbl_In.ASCImgTimeStart_s(GroupID_AllVidData_WBCConcStability_Tbl_In.ASCImgTimeStart_s>eps),[],1)]);
title({strcat('ASCRatioDiffLim:',32,num2str(ASCRatioDiffLim)),...
    strcat('# GroupID:',32,num2str(height(bpData)),32,'(out of',32,num2str(height(GroupID_AllVidData_WBCConcStability_Tbl_In)),')'),...
    strcat('Median:',32,num2str(median(bpData))),...
    strcat('25-Percentile:',32,num2str(quantile(bpData,0.25))),...
    strcat('75-Percentile:',32,num2str(quantile(bpData,0.75))),...
    strcat('Upper whisker:',32,num2str(quantile(bpData,0.75)+1.5*(quantile(bpData,0.75)-quantile(bpData,0.25)))),...
    strcat('Lower whisker:',32,num2str(quantile(bpData,0.25)-1.5*(quantile(bpData,0.75)-quantile(bpData,0.25)))),...
    strcat('Mean:',32,num2str(mean(bpData))),...
    strcat('Stdev:',32,num2str(std(bpData))),...
    strcat('FlowVolBinSize:',32,num2str(min(ASC_ReportIncVol_Step_In)),'-',num2str(max(ASC_ReportIncVol_Step_In)),'nL'),...
    },...
    'FontSize',7,'FontWeight','normal','Interpreter','none');


ax2 = nexttile(tlo); 
hold(ax2, 'on');
grid on;
bpData = GroupID_AllVidData_WBCConcStability_Tbl_In.ASCCumFlowVolStart_nL(GroupID_AllVidData_WBCConcStability_Tbl_In.ASCCumFlowVolStart_nL>eps);
boxplot(ax2, bpData);
xlabel('All_GroupID','interpreter','none');
xlim([0.5 1.5]);
if UI_In.ASCOption==1
    ylabel('CumFlowVolStart_nL (IncrWBCConc)','interpreter','none');
elseif UI_In.ASCOption==2
    ylabel('CumFlowVolStart_nL (CumWBCConc)','interpreter','none')
end
ylim([0 max(GroupID_AllVidData_WBCConcStability_Tbl_In.ASCCumFlowVolStart_nL(GroupID_AllVidData_WBCConcStability_Tbl_In.ASCCumFlowVolStart_nL>eps),[],1)]);
title({strcat('ASCRatioDiffLim:',32,num2str(ASCRatioDiffLim)),...    
    strcat('# GroupID:',32,num2str(height(bpData)),32,'(out of',32,num2str(height(GroupID_AllVidData_WBCConcStability_Tbl_In)),')'),...
    strcat('Median:',32,num2str(median(bpData))),...
    strcat('25-Percentile:',32,num2str(quantile(bpData,0.25))),...
    strcat('75-Percentile:',32,num2str(quantile(bpData,0.75))),...
    strcat('Upper whisker:',32,num2str(quantile(bpData,0.75)+1.5*(quantile(bpData,0.75)-quantile(bpData,0.25)))),...
    strcat('Lower whisker:',32,num2str(quantile(bpData,0.25)-1.5*(quantile(bpData,0.75)-quantile(bpData,0.25)))),...
    strcat('Mean:',32,num2str(mean(bpData))),...
    strcat('Stdev:',32,num2str(std(bpData))),...
    strcat('FlowVolBinSize:',32,num2str(min(ASC_ReportIncVol_Step_In)),'-',num2str(max(ASC_ReportIncVol_Step_In)),'nL'),...
    },...
    'fontsize',7,'fontweight','normal','interpreter','none');

% = Plot Output

Plot_ASCVolBinStats_Out = export_fig('-r150');

close all;
clearvars ASC ASCRatioDiffLim bpData FigHandle ax* tlo VesselDiameterHistEdge VesselDiameterHistBinCtr grp_VesselDiameterHistBin grpID_VesselDiameterHistBin;
