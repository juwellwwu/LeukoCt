function [Plot_ASCVolBinConc_gCell_Out] = fPlotASCVolBinConc(GroupID_AllVidData_gCell_In, UI_In)

if UI_In.ASCOption == 1 % Incremental Flow Vol Bin
    ASCRatioDiffLim_In = UI_In.ASC.IncrVolBin.ASCRatioDiffLim;
elseif UI_In.ASCOption == 2
    ASCRatioDiffLim_In = UI_In.ASC.CumVolBin.ASCRatioDiffLim;
end

%% Plot

Plot_ASCVolBinConc_gCell_Out = cell(height(GroupID_AllVidData_gCell_In),1);

for ugCt = 1:height(GroupID_AllVidData_gCell_In)

    BestVolStep_ASCParamRow_In = find(GroupID_AllVidData_gCell_In{ugCt,5}(:,10));

    if isempty(BestVolStep_ASCParamRow_In) | all(GroupID_AllVidData_gCell_In{ugCt,5}(:,10)<eps) 
        Plot_ASCVolBinConc_gCell_Out{ugCt,1} = [];
    else
        BestVolStep_RData_In = GroupID_AllVidData_gCell_In{ugCt,4}{1,1,BestVolStep_ASCParamRow_In};
        BestVolStep_CumVolStart_In = GroupID_AllVidData_gCell_In{ugCt,8};
        BestVolStep_In = GroupID_AllVidData_gCell_In{ugCt,9};
        BestVolStep_BinIdxStart_In = GroupID_AllVidData_gCell_In{ugCt,10};
        BestVolStep_BinIdxEnd_In = GroupID_AllVidData_gCell_In{ugCt,11};
        BestVolStep_MaxRatioDiff_In = GroupID_AllVidData_gCell_In{ugCt,12};

        plData = GroupID_AllVidData_gCell_In{ugCt,2};

        % = Plot

        FigHandle = figure('Position',[10 50 1500 800],'visible','off');

        tlo = tiledlayout(4,1,'TileSpacing','Compact','Padding','tight');
        cmap1 = colormap(lines(numel(unique(plData(:,1,1)))));
        cmap2 = cmap1.*0.67;

        % = (Incremental) Blood flow (nL) at tMidPt_FrCtr (sec), referenced to start of video 1 (variation among Orient Bands)
        ax = nexttile(tlo);
        hold(ax, 'on');
        grid on;
        gscatter(plData(:,5),plData(:,6),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off'); % Color by video
        gscatter(plData(:,5),plData(:,6)-plData(:,7),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off');
        gscatter(plData(:,5),plData(:,6)+plData(:,7),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off');
        plot(plData(:,5),plData(:,6),'LineWidth',1,'Color',[0 0 0]);
        xlabel('tMidPt_FrCtr (s)','Interpreter','none');
        ylabel('IncrFlowVol_OB (nL)','Interpreter','none');
        xlim([0 10*ceil(plData(end,3)/10)]);
        xticks(0:10:10*ceil(plData(end,3)/10));
        set(gca,'TickDir','out');

        title({strcat('GroupID:',32,GroupID_AllVidData_gCell_In{ugCt,1})}, ...
            'FontSize',10,'FontWeight','normal','Interpreter','none');

        % = Cumulative WBC conc (#/nL = K/uL) at tMidPt_FrCtr (sec), referenced to start of video 1 (variation among Orient Bands)
        ax2 = nexttile(tlo);
        hold(ax2, 'on');
        grid on;
        gscatter(plData(:,3),plData(:,37),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off'); % Color by video
        gscatter(plData(:,3),plData(:,38),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off');
        gscatter(plData(:,3),plData(:,39),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off');
        plot(plData(:,3),plData(:,37),'LineWidth',1,'Color',[0 0 0]);
        xlabel('tMidPt_FrEnd (s)','Interpreter','none');
        ylabel('CumWBCConc_OB (K/uL)','Interpreter','none');
        xlim([0 10*ceil(plData(end,3)/10)]);
        xticks(0:10:10*ceil(plData(end,3)/10));
        ylim([floor(min(plData(:,38),[],1)) ceil(max(plData(:,39),[],1))]);
        set(gca,'TickDir','out');

        % = Cumulative Blood flow (nL) vs Cumulative WBC conc (#/nL = K/uL) (variation among Orient Bands)
        ax3 = nexttile(tlo); 
        hold(ax3, 'on');
        grid on;
        gscatter(plData(:,31,1),plData(:,37,1),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off'); 
        gscatter(plData(:,31,1),plData(:,38,1),plData(:,1),cmap1,'.',repelem(2,numel(unique(plData(:,1)))),'doleg','off'); 
        gscatter(plData(:,31,1),plData(:,39,1),plData(:,1),cmap1,'.',repelem(2,numel(unique(plData(:,1)))),'doleg','off'); 
        plot(plData(:,31,1),plData(:,37,1),'LineWidth',1,'Color',[0 0 0]);
        xline(BestVolStep_RData_In(:,2),":",'LineWidth',1);
        xline(BestVolStep_CumVolStart_In,"--",'LineWidth',2); 
        xlabel('CumFlowVol_OB (nL)','Interpreter','none');
        ylabel('CumWBCConc_OB (K/uL)','Interpreter','none');
        xlim([0 5*ceil(plData(end,31,1)/5)]);
        xticks(0:5:5*ceil(plData(end,31,1)/5));
        ylim([floor(min(plData(:,38,1),[],1)) ceil(max(plData(:,39,1),[],1))]);
        set(gca,'TickDir','out');
        title({strcat('Cumulative Flow Vol (nL) Start for Stable WBC Conc (', num2str(BestVolStep_In), 'nL Flow Vol Bins):',32,num2str(BestVolStep_CumVolStart_In,'%0.2f'),'nL')}, ...
            'FontSize',10,'FontWeight','normal','Interpreter','none');

        % = Line Plot describing concentration
        ax4 = nexttile(tlo);
        hold(ax4, 'on');
        grid on;
        cmap2 = colormap(jet(2));
        x_plot = sort(vertcat(BestVolStep_RData_In(:,2),BestVolStep_RData_In(:,3)),'ascend');
        y_plot = repelem(BestVolStep_RData_In(:,7),2);
        for lnCt = 1:height(BestVolStep_RData_In)
            if ismember(lnCt,(BestVolStep_BinIdxStart_In:BestVolStep_BinIdxEnd_In)) % Stable WBC Bin
                lnColor = cmap2(2,:);
            else
                lnColor = cmap2(1,:);
            end
            lnHandle = plot(BestVolStep_RData_In(lnCt,2):UI_In.CumVolIntp_Step:BestVolStep_RData_In(lnCt,3),... % x plot data
                repelem(BestVolStep_RData_In(lnCt,7),numel(BestVolStep_RData_In(lnCt,2):UI_In.CumVolIntp_Step:BestVolStep_RData_In(lnCt,3))),... % y plot data
                'LineWidth',5,'Color',lnColor);
            label(lnHandle,num2str(BestVolStep_RData_In(lnCt,7)),'location','top'); 
            label(lnHandle,num2str(BestVolStep_RData_In(lnCt,9)),'location','bottom'); 
        end
        xline(BestVolStep_RData_In(:,2),":",'LineWidth',1); 
        xline(BestVolStep_CumVolStart_In,"--",'LineWidth',2); 
        xlim([0 5*ceil(plData(end,31)/5)]);
        xticks(0:5:5*ceil(plData(end,31)/5));
        xlabel('CumFlowVol_OB (nL)','Interpreter','none');
        if UI_In.ASCOption == 1
            ylabel('Bin Incr WBC Conc (K/uL)','Interpreter','none');
        elseif UI_In.ASCOption == 2
            ylabel('Bin Cum WBC Conc (K/uL)','Interpreter','none');
        end
        set(gca,'TickDir','out');
        title({strcat('Max. WBC Conc Ratio-Difference-from-1 among Stable WBC Conc Flow Vol Bins:', num2str(BestVolStep_MaxRatioDiff_In)), ...
            strcat('User-defined, Max. Allowed WBC Conc Ratio-Difference-from-1:',32,num2str(ASCRatioDiffLim_In))}, ...
            'FontSize',10,'FontWeight','normal','Interpreter','none');

        % = Load into Output
        Plot_ASCVolBinConc_gCell_Out{ugCt,1} = export_fig('-r150');

        close all;
        clearvars FigHandle tlo ax* cmap* x_plot y_plot lnCt lnColor;

    end

end