function Fig_Out = fPlotFrIntensityProfile(FrIntensitySum_In,Fr_AutoCellCt_Struc_In,Fr_ManualCellCt_In,tFr_Flag_fStruct_In,text_yIncrement_In,yScale_In)

FrIntensitySum_In = FrIntensitySum_In.*yScale_In;
NumBlk_In = width(FrIntensitySum_In);

if ~isempty(Fr_AutoCellCt_Struc_In) & min(Fr_AutoCellCt_Struc_In.CellCt,[],'all')>-eps
    AutoCellCt_In = Fr_AutoCellCt_Struc_In.CellCt;
    AutoCellCt_pp_In = Fr_AutoCellCt_Struc_In.pkParam; 
    AutoCellCt_Rank_In = Fr_AutoCellCt_Struc_In.BrightRank; 
    AutoCellCt_Optmz_In = Fr_AutoCellCt_Struc_In.Optimized; 
end

if ~isempty(tFr_Flag_fStruct_In) 

    tFr_Flag_fStruct_In = orderfields(tFr_Flag_fStruct_In); % order fields
    tFr_Flag_fStruct_Fields_In = fieldnames(tFr_Flag_fStruct_In);
    FlagID = cellfun(@isempty,strfind(tFr_Flag_fStruct_Fields_In, 'Flag')); 
    tFr_Flag_fStruct_Fields_In(FlagID) = [];

    tFr_Flag_PlotData_fCell = cell(height(fieldnames(tFr_Flag_fStruct_In)),1);
    tFr_Flag_PlotColor_fCell = cell(height(fieldnames(tFr_Flag_fStruct_In)),1);

    tFr_Flag_Ploty = ... % height of Flag lines
        linspace(min(FrIntensitySum_In,[],'all')-0.05*range(FrIntensitySum_In,'all'),min(FrIntensitySum_In,[],'all')-0.20*range(FrIntensitySum_In,'all'),...
        height(tFr_Flag_fStruct_Fields_In));

    for fCt = 1:height(tFr_Flag_fStruct_Fields_In) 
        tFr_Flag_Plot_fCell{fCt,1} = getfield(tFr_Flag_fStruct_In,tFr_Flag_fStruct_Fields_In{fCt});
        tFr_Flag_PlotData_fCell{fCt,1} = single(tFr_Flag_Plot_fCell{fCt,1}.Data); 
        tFr_Flag_PlotData_fCell{fCt,1}(tFr_Flag_PlotData_fCell{fCt,1}<eps) = NaN; 
        tFr_Flag_PlotData_fCell{fCt,1} = tFr_Flag_PlotData_fCell{fCt,1}.*tFr_Flag_Ploty(fCt);
        if width(tFr_Flag_PlotData_fCell{fCt,1})<NumBlk_In
            tFr_Flag_PlotData_fCell{fCt,1} = repmat(tFr_Flag_PlotData_fCell{fCt,1},1,NumBlk_In);
        end

        tFr_Flag_PlotColor_fCell{fCt,1} = tFr_Flag_Plot_fCell{fCt,1}.Color; 
    end

end

cmap = colormap(summer(round(NumBlk_In*1.5))); % Vessel blocks: Dark green to grass green

FigHandle = figure('Position',[100 50 1600 NumBlk_In*300],'visible','on');

if ~isempty(Fr_ManualCellCt_In)
    tlo = tiledlayout(NumBlk_In+1,1,'TileSpacing','Compact','Padding','tight'); 
    ax_Cell = cell(NumBlk_In+1,1);
else
    tlo = tiledlayout(NumBlk_In,1,'TileSpacing','Compact','Padding','tight');
    ax_Cell = cell(NumBlk_In,1);
end

for BlkCt = 1:NumBlk_In

    ax_Cell{BlkCt,1} = nexttile(tlo);
    hold(ax_Cell{BlkCt,1}, 'on');
    grid on;

    if FrIntensitySum_In(1,BlkCt)>(-9999) 

        if ~isempty(Fr_AutoCellCt_Struc_In) & min(Fr_AutoCellCt_Struc_In.CellCt,[],'all')>-eps

            AutoCellCt_PeakHeight = AutoCellCt_pp_In{1,BlkCt}(:,1).*yScale_In; 
            AutoCellCt_PeakLoc = AutoCellCt_pp_In{1,BlkCt}(:,2); 

        else

            AutoCellCt_PeakHeight = [];
            AutoCellCt_PeakLoc = [];
            AutoCellCt_In = -9999;

        end

        % = Line Plot
        plot((1:1:height(FrIntensitySum_In)),FrIntensitySum_In(:,BlkCt),'Color',cmap(BlkCt,:),'Linewidth',1.0);

        if ~isempty(Fr_ManualCellCt_In)
            xline(Fr_ManualCellCt_In,':','Color','r','Linewidth',0.50);
        end

        if ~isempty(tFr_Flag_fStruct_In)
            for fCt = 1:height(tFr_Flag_fStruct_Fields_In) % # different flags
                hold on;
                plot((1:1:height(FrIntensitySum_In)),tFr_Flag_PlotData_fCell{fCt,1}(:,BlkCt),'Color',tFr_Flag_PlotColor_fCell{fCt,1},'Linewidth',1.50)
            end
        end

        if ~isempty(AutoCellCt_PeakLoc)

            AutoCellCt_PeakLocMark = repelem('*',height(AutoCellCt_PeakLoc))';
            text(AutoCellCt_PeakLoc,AutoCellCt_PeakHeight+0.005,AutoCellCt_PeakLocMark,...
                'HorizontalAlignment','center','FontSize',3,'Color',[0 0.5 0]);

        end

        % % Set Axis Limits
        if ~isempty(tFr_Flag_fStruct_In)
            ymin =min(min(FrIntensitySum_In(:,BlkCt),[],'all')-0.02*range(FrIntensitySum_In(:,BlkCt),'all'),min(tFr_Flag_Ploty,[],"all")-0.02*range(FrIntensitySum_In(:,BlkCt),'all'));
        else
            ymin = min(FrIntensitySum_In(:,BlkCt),[],'all')-0.02*range(FrIntensitySum_In(:,BlkCt),'all');
        end
        if ~isempty(AutoCellCt_PeakLoc)
            ymax = max(AutoCellCt_PeakHeight,[],'all')+0.01; 
        else
            ymax = max(FrIntensitySum_In(:,BlkCt),[],'all')+0.02*range(FrIntensitySum_In(:,BlkCt),'all');
        end
        ylim([ymin ymax]);

        xlim([0 height(FrIntensitySum_In)]);
        xticks(0:100:100*floor(0.01*height(FrIntensitySum_In)));

        xlabel('Frame');
        ylabel(strcat('FrIntensitySum'));
        fontsize(gca,10,"pixels");

        if ~isempty(Fr_ManualCellCt_In) & (~isempty(Fr_AutoCellCt_Struc_In) & min(Fr_AutoCellCt_Struc_In.CellCt,[],'all')>-eps)
            title({strcat('AutoCellCt:',32,num2str(AutoCellCt_In(1,BlkCt)),'; ManualCellCt:',32,num2str(numel(Fr_ManualCellCt_In)),'; ManualCellCt (no repeat):',32,num2str(numel(unique(Fr_ManualCellCt_In))));...
                strcat('BrightRank:',32,num2str(AutoCellCt_Rank_In(1,BlkCt)),',',32,'Optimized:',32,num2str(AutoCellCt_Optmz_In(1,BlkCt)))},...
                'FontSize',8,'FontWeight','normal');
        elseif ~isempty(Fr_ManualCellCt_In) & (isempty(Fr_AutoCellCt_Struc_In) | min(Fr_AutoCellCt_Struc_In.CellCt,[],'all')<0)
            title({strcat('ManualCellCt:',32,num2str(numel(Fr_ManualCellCt_In)),'; ManualCellCt (no repeat):',32,num2str(numel(unique(Fr_ManualCellCt_In))));...
                strcat('BrightRank: NA')},...
                'FontSize',8,'FontWeight','normal');
        elseif isempty(Fr_ManualCellCt_In) & (~isempty(Fr_AutoCellCt_Struc_In) & min(Fr_AutoCellCt_Struc_In.CellCt,[],'all')>-eps)
            title({strcat('AutoCellCt:',32,num2str(AutoCellCt_In(1,BlkCt)));...
                strcat('BrightRank:',32,num2str(AutoCellCt_Rank_In(1,BlkCt)),',',32,'Optimized:',32,num2str(AutoCellCt_Optmz_In(1,BlkCt)))},...
                'FontSize',8,'FontWeight','normal');
        end

    else 

        xlim([0 height(FrIntensitySum_In)]);
        xticks(0:100:100*floor(0.01*height(FrIntensitySum_In)));
        xlabel('Frame');
        ylabel(strcat('FrIntensitySum'));
        fontsize(gca,10,"pixels");
        title({strcat('(No Cell Count: y-data not given in Fr_ManualCellCt_In)')},...
            'FontSize',8,'FontWeight','normal','Interpreter','none');

    end

end

% If manual ct information available, add plot with info
if ~isempty(Fr_ManualCellCt_In)

    ax_Cell{NumBlk_In+1,1} = nexttile(tlo);
    hold(ax_Cell{NumBlk_In+1,1}, 'on');
    grid on;

    xline(Fr_ManualCellCt_In,':','Color','r','Linewidth',0.50);

    if ~isempty(Fr_ManualCellCt_In)
        text(Fr_ManualCellCt_In(1:8:end,1),repelem(0+0.50,numel(Fr_ManualCellCt_In(1:8:end,1))),num2str(Fr_ManualCellCt_In(1:8:end,1)),...
            'HorizontalAlignment','center','FontSize',6,'Color',[0.5 0 0]);
        text(Fr_ManualCellCt_In(2:8:end,1),repelem(0+0.5+1*text_yIncrement_In,numel(Fr_ManualCellCt_In(2:8:end,1))),num2str(Fr_ManualCellCt_In(2:8:end,1)),...
            'HorizontalAlignment','center','FontSize',6,'Color',[0.5 0 0]);
        text(Fr_ManualCellCt_In(3:8:end,1),repelem(0+0.5+2*text_yIncrement_In,numel(Fr_ManualCellCt_In(3:8:end,1))),num2str(Fr_ManualCellCt_In(3:8:end,1)),...
            'HorizontalAlignment','center','FontSize',6,'Color',[0.5 0 0]);
        text(Fr_ManualCellCt_In(4:8:end,1),repelem(0+0.5+3*text_yIncrement_In,numel(Fr_ManualCellCt_In(4:8:end,1))),num2str(Fr_ManualCellCt_In(4:8:end,1)),...
            'HorizontalAlignment','center','FontSize',6,'Color',[0.5 0 0]);
        text(Fr_ManualCellCt_In(5:8:end,1),repelem(0+0.5+4*text_yIncrement_In,numel(Fr_ManualCellCt_In(5:8:end,1))),num2str(Fr_ManualCellCt_In(5:8:end,1)),...
            'HorizontalAlignment','center','FontSize',6,'Color',[0.5 0 0]);
        text(Fr_ManualCellCt_In(6:8:end,1),repelem(0+0.5+5*text_yIncrement_In,numel(Fr_ManualCellCt_In(6:8:end,1))),num2str(Fr_ManualCellCt_In(6:8:end,1)),...
            'HorizontalAlignment','center','FontSize',6,'Color',[0.5 0 0]);
        text(Fr_ManualCellCt_In(7:8:end,1),repelem(0+0.5+6*text_yIncrement_In,numel(Fr_ManualCellCt_In(7:8:end,1))),num2str(Fr_ManualCellCt_In(7:8:end,1)),...
            'HorizontalAlignment','center','FontSize',6,'Color',[0.5 0 0]);
        text(Fr_ManualCellCt_In(8:8:end,1),repelem(0+0.5+7*text_yIncrement_In,numel(Fr_ManualCellCt_In(8:8:end,1))),num2str(Fr_ManualCellCt_In(8:8:end,1)),...
            'HorizontalAlignment','center','FontSize',6,'Color',[0.5 0 0]);
    end

    xlim([0 height(FrIntensitySum_In)]);
    xticks(0:100:100*floor(0.01*height(FrIntensitySum_In)));

    ylim([0 30]);

end

Fig_Out = export_fig(FigHandle,'-a1','-r150');

end
