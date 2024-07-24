function [Fig_VidCombotMidPtData_Cell_Out] = fPlotVidCombotMidPtData(GroupID_AllVidData_gCell_In)

Fig_VidCombotMidPtData_Cell_Out = cell(height(GroupID_AllVidData_gCell_In),1);

fprintf(strcat('Plotting all-video-combined, flag-cut tMidPt flow volume and WBC conc data...\n'));

for ugCt = 1:height(GroupID_AllVidData_gCell_In) % unique GroupID

    plData = GroupID_AllVidData_gCell_In{ugCt,2}; % plot data

    if isempty(plData)

        fprintf(strcat('GroupID',32,num2str(ugCt),'/',num2str(height(GroupID_AllVidData_gCell_In)),': No unflagged tMidPt','...\n'));

        Fig_VidCombotMidPtData_Cell_Out{ugCt,1} = [];
    
    else
    
        fprintf(strcat('GroupID',32,num2str(ugCt),'/',num2str(height(GroupID_AllVidData_gCell_In)),32,'...\n'));

        GroupIDStr = GroupID_AllVidData_gCell_In{ugCt,1};

        FigHandle = figure('Position',[20 100 1500 400],'visible','off');

        tlo = tiledlayout(3,1,'TileSpacing','Compact','Padding','tight');
        cmap1 = colormap(lines(numel(unique(plData(:,1))))); % # videos
        cmap2 = cmap1.*0.67;

        % = Incremental Blood Volume (nL) at tMidPt_FrCtr (sec), referenced to start of video 1 (variation among Orient Bands OB)
        ax = nexttile(tlo);
        hold(ax, 'on');
        grid on;
        gscatter(plData(:,5),plData(:,6),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off'); 
        gscatter(plData(:,5),plData(:,6)-plData(:,7),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off');
        gscatter(plData(:,5),plData(:,6)+plData(:,7),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off');
        plot(plData(:,5),plData(:,6),'LineWidth',1,'Color',[0 0 0]);
        xlabel('tMidPt_FrCtr (s)','Interpreter','none');
        ylabel('IncrFlowVol_OB (nL)','Interpreter','none'); 
        xlim([0 10*ceil(plData(end,3)/10)]);
        xticks(0:10:10*ceil(plData(end,3)/10));
        set(gca,'TickDir','out');

        title({strcat('GroupID:',32,GroupIDStr)}, ...
            'FontSize',10,'FontWeight','normal','Interpreter','none');

        % = Incremental WBC Conc (K/uL) at tMidPt_FrCtr (sec), referenced to start of video 1 (variation among Orient Bands OB)      
        ax3 = nexttile(tlo);
        hold(ax3, 'on');
        grid on;
        gscatter(plData(:,5),plData(:,12),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off'); % Color by video
        gscatter(plData(:,5),plData(:,13),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off');
        gscatter(plData(:,5),plData(:,14),plData(:,1),cmap1,'.',repelem(6,numel(unique(plData(:,1)))),'doleg','off');
        plot(plData(:,5),plData(:,12),'LineWidth',1,'Color',[0 0 0]);
        xlabel('tMidPt_FrCtr (s)','Interpreter','none');
        ylabel('IncrWBCConc_OB (K/uL)','Interpreter','none');
        xlim([0 10*ceil(plData(end,3)/10)]);
        xticks(0:10:10*ceil(plData(end,3)/10));
        set(gca,'TickDir','out');

        % = Cumulative WBC conc (#/nL = K/uL) at tMidPt_FrCtr (sec), referenced to start of video 1 (variation among Orient Bands OB)
        ax5 = nexttile(tlo);
        hold(ax5, 'on');
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

        fontsize(gcf,scale=0.7);

        Fig_VidCombotMidPtData_Cell_Out{ugCt,1} = export_fig('-r150');

        close all;

        clearvars cmap* tlo ax* FigHandle;

    end

end

clearvars ugCt plData GroupIDStr;
% % % clearvars GrpID_Idx GrpID_uq ugCt PipeOut_FolderFileList;
