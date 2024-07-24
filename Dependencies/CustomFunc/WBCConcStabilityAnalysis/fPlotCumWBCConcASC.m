function [Plot_GroupID_CumTorVvsCumConc_lpmCell_Out, PlotLabel_GroupID_CumTorVvsCumConc_lpmCell_Out] = ...
    fPlotCumWBCConcASC(GroupID_CumTorVvsCumConc_PlotData_torvlCell_In,UI_In,TimeorVolOption_In)

%% Designate all Cumulative Imaging Time or Flow Volume Limits for plotting 

if TimeorVolOption_In == 1 
    xStr = 'ImgTime_s';
    CumTorVLim_All = UI_In.CumTimeLim_All;
elseif TimeorVolOption_In ==2 
    xStr = 'FlowVol_nL';
    CumTorVLim_All = UI_In.CumVolLim_All;
end
if ~iscolumn(CumTorVLim_All)
    CumTorVLim_All = CumTorVLim_All';
end


%% Designate Normalization Option

NormOption = UI_In.CumTimeVolvsCumConcNormOption;
NormStr = strcat('_Norm',num2str(NormOption));


%% Plot 

Plot_GroupID_CumTorVvsCumConc_lpmCell_Out = cell(numel(CumTorVLim_All),2,2); 
PlotLabel_GroupID_CumTorVvsCumConc_lpmCell_Out = cell(numel(CumTorVLim_All),2,2); 

for mCt = 1:2 

    if mCt == 1
        MethodStr = '_OB';
    elseif mCt == 2
        MethodStr = '_SK';
    end
 
    for torvlCt = 1:numel(CumTorVLim_All) 

        CumTorVLim = CumTorVLim_All(torvlCt); 
        GroupID_CumTorVvsCumConc_Trim = GroupID_CumTorVvsCumConc_PlotData_torvlCell_In{mCt,1}{torvlCt,7};
        GroupID_Label = GroupID_CumTorVvsCumConc_PlotData_torvlCell_In{mCt,1}{torvlCt,8}.GroupID;
        WinVesselDiameter_um = GroupID_CumTorVvsCumConc_PlotData_torvlCell_In{mCt,1}{torvlCt,8}.WinVesselDiameter_um;

        % = Plot
        Fig_CumTorVvsCumWBCConc_pCell = cell(2,1); 

        for pCt = 1:2 

            if pCt == 1
                pStr = '_wLoHi';
            elseif pCt == 2
                pStr = '_NoLoHi';
            end

            [~,~,VesselDiameterHistBin] = histcounts(WinVesselDiameter_um,UI_In.VesselDiameterHistEdge_um);
            cmap = colormap(jet(numel(UI_In.VesselDiameterHistEdge_um)-1)); % Rainbow; Small diameter = Blue; Large diameter = Red
            cmap = vertcat(cmap,[0.8 0.8 0.8]); 
            VesselDiameterHistBin(VesselDiameterHistBin==0) = numel(UI_In.VesselDiameterHistEdge_um); 
            cmap = cmap(VesselDiameterHistBin,:);
            
            close all;

            FigHandle = figure('Position',[10 50 1500 800],'visible','off');
            tlo = tiledlayout(1,1,'TileSpacing','Compact','Padding','tight');

            ax = nexttile(tlo);
            hold(ax, 'on');
            grid on;

            plot(GroupID_CumTorVvsCumConc_Trim(:,1,1),GroupID_CumTorVvsCumConc_Trim(:,2:end,1),'LineWidth',1.5); 
            if pCt == 1
                plot(GroupID_CumTorVvsCumConc_Trim(:,1,2),GroupID_CumTorVvsCumConc_Trim(:,2:end,2),'--','LineWidth',0.5); 
                plot(GroupID_CumTorVvsCumConc_Trim(:,1,3),GroupID_CumTorVvsCumConc_Trim(:,2:end,3),'--','LineWidth',0.5); 
            end
            colororder(cmap);
            legend(GroupID_Label,'Orientation','vertical','Location','eastoutside');
            xlabel(strcat('Cum',xStr),'Interpreter','none');
            ylabel('Cum WBC Conc (K/uL)','Interpreter','none');
            xlim([0 CumTorVLim]);
            xticks(0:10:CumTorVLim);
            ylim([0.5*floor(min(GroupID_CumTorVvsCumConc_Trim(:,2:end,2),[],'all')/0.5) 0.5*ceil(max(GroupID_CumTorVvsCumConc_Trim(ceil(0.01*height(GroupID_CumTorVvsCumConc_Trim)):end,2:end,3),[],'all')/0.5)]);
            yticks([0.5*floor(min(GroupID_CumTorVvsCumConc_Trim(:,2:end,2),[],'all')/0.5):0.5:0.5*ceil(max(GroupID_CumTorVvsCumConc_Trim(ceil(0.01*height(GroupID_CumTorVvsCumConc_Trim)):end,2:end,3),[],'all')/0.5)]);
     
            set(gca,'TickDir','out');

            % = Save Plot
            Plot_GroupID_CumTorVvsCumConc_lpmCell_Out{torvlCt,pCt,mCt} = export_fig('-r150');
            PlotLabel_GroupID_CumTorVvsCumConc_lpmCell_Out{torvlCt,pCt,mCt} = strcat(NormStr,MethodStr,pStr,'_0-',num2str(CumTorVLim));

            close all;
            clearvars FigHandle ax cmap tlo;

        end
        
    end

  end


close all;
clearvars FigHandle ax NormOption NormStr tlCt CumTorVLim_All ugCt2 pCt mCt MethodStr xStr;
clearvars GroupID_CumTorVvsCumConc GroupID_CumTorVvsCumConc_PlotData_torvlCell;

