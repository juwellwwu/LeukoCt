function [GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_Out, ...
    Plot_GroupID_CumTorVvsCumConc_lpmCell_Out, PlotLabel_GroupID_CumTorVvsCumConc_lpmCell_Out, ...
    Tbl_GroupID_CumTorVvsCumConc_lmCell_Out, TblLabel_GroupID_CumTorVvsCumConc_lmCell_Out] = ...
    fPlotCumWBCConc(GroupID_CumTorVvsCumConc_mCell_In,GroupID_AllVidData_gCell_In,UI_In,TimeorVolOption_In)

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

% = Cell for holding plot data
GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_Out = cell(2,1);

% = Cell for holding plot image
Plot_GroupID_CumTorVvsCumConc_lpmCell_Out = cell(numel(CumTorVLim_All),2,2); 
PlotLabel_GroupID_CumTorVvsCumConc_lpmCell_Out = cell(numel(CumTorVLim_All),2,2); 

% = Cell for holding Excel CumTorV vs CumConc Table
Tbl_GroupID_CumTorVvsCumConc_lmCell_Out = cell(numel(CumTorVLim_All),1,2); 
TblLabel_GroupID_CumTorVvsCumConc_lmCell_Out = cell(numel(CumTorVLim_All),1,2);  

for mCt = 1:2 

    if mCt == 1
        MethodStr = '_OB';
    elseif mCt == 2
        MethodStr = '_SK';
    end
    GroupID_CumTorVvsCumConc = GroupID_CumTorVvsCumConc_mCell_In{mCt,1};

    % = Prepare "GroupID_CumTorVvsCumConc_PlotData_torvlCell " cell to hold plot data
    % "TorV": T = Cumulative Imaging Time (s); V = Cumulative Flow Volume (nL)
    % # Rows = # CumTorVLim
    % Col 1) CumTorVLim 
    % Col 2) Normalization Method (0-4)
    % Col 3) GroupID_CumTorVvsCumConc Data for GroupID w/ ImgTime -OR- Flow Vol >=
    % CumTorVLimit (= GroupID_CumTorVvsCumConc_Trim)
    % Col 4) Window Vessel Diameter (um), REFCellConc (K/uL), for GroupID w/ ImgTime -OR- Flow Vol >=
    % CumTorVLimit (= WinDiameterCBCTbl)
    % Col 5) ASC WBC Conc Stability Option: 1 if use Incremental Vol Bin, 2 if use
    % Cumulative Vol Bin (= UI.ASCOption)
    % Col 6) ASC WBC Conc Stability: Max allowed WBC Conc
    % Ratio-Difference-from-1 (= UI.ASCRatioDiffLim_IncrVolBin or = UI.ASCRatioDiffLim_CumVolBin)
    % Col 7) GroupID_CumTorVvsCumConc Data for GroupID w/ ImgTime -OR- Flow Vol >=
    % CumTorVLimit & Pass ASC Stabillty
    % Col 8) Window Vessel Diameter (um), REFCellConc (K/uL), for GroupID w/ ImgTime -OR- Flow Vol >=
    % CumTorVLimit (= WinDiameterCBCTbl) & Pass ASC Stabillty

    GroupID_CumTorVvsCumConc_PlotData_torvlCell = cell(numel(CumTorVLim_All),8); % torvl = Time or Vol Limit

    for torvlCt = 1:numel(CumTorVLim_All) % VolLim Ct

        % = Restrict analysis to GroupID w/ Cum Time or Flow Vol >= CumTorVLim
        CumTorVLim = CumTorVLim_All(torvlCt); 
        GroupID_CumTorVvsCumConc_Trim = GroupID_CumTorVvsCumConc;
        if NormOption==0 | NormOption==1 | NormOption==2
            FlagID = isnan(GroupID_CumTorVvsCumConc_Trim(GroupID_CumTorVvsCumConc_Trim(:,1,1)==CumTorVLim+CumTorVIntp_Step_In,:,1)); 
            GroupID_CumTorVvsCumConc_Trim(:,FlagID,:) = []; 
        elseif NormOption==3 | NormOption==4 
            FlagID = isnan(GroupID_CumTorVvsCumConc_Trim(1,:,1)); 
            GroupID_CumTorVvsCumConc_Trim(:,FlagID,:) = []; 
        end
        GroupID_CumTorVvsCumConc_Trim(GroupID_CumTorVvsCumConc_Trim(:,1)>CumTorVLim,:,:) = [];

        % = Option to Normalize Data
        if NormOption == 1
            GroupID_CumTorVvsCumConc_Trim(:,2:end,:) = ... 
                GroupID_CumTorVvsCumConc_Trim(:,2:end,:)./GroupID_CumTorVvsCumConc_Trim(end,2:end,1);
        elseif NormOption == 2
            GroupID_CumTorVvsCumConc_Trim(:,2:end,:) = ... 
                GroupID_CumTorVvsCumConc_Trim(:,2:end,:)./GroupID_CumTorVvsCumConc_Trim(end,2:end,1);
            GroupID_CumTorVvsCumConc_Trim(:,2:end,:) = ... 
                abs(GroupID_CumTorVvsCumConc_Trim(:,2:end,:)-1);
        elseif NormOption == 3 
            CumConc_LastNonNaNIdx = ...
                arrayfun(@(x) find(~isnan(GroupID_CumTorVvsCumConc_Trim(:,x,1)),1,'last'),1:width(GroupID_CumTorVvsCumConc_Trim));
        elseif NormOption == 4 
            CumConc_LastNonNaNIdx = ... 
                arrayfun(@(x) find(~isnan(GroupID_CumTorVvsCumConc_Trim(:,x,1)),1,'last'),1:width(GroupID_CumTorVvsCumConc_Trim));
            for ugCt2 = 2:width(GroupID_CumTorVvsCumConc_Trim) 
                GroupID_CumTorVvsCumConc_Trim(:,ugCt2,:) = ... 
                    GroupID_CumTorVvsCumConc_Trim(:,ugCt2,:)./GroupID_CumTorVvsCumConc_Trim(CumConc_LastNonNaNIdx(ugCt2),ugCt2,1);
            end
        end

        % = Plot
        Fig_CumTorVvsCumWBCConc_pCell = cell(2,1); % 2 plot choices: one with WBC Conc Lo and Hi lines, one without

        for pCt = 1:2 % one w/ Lo and Hi lines, one w/out

            if pCt == 1
                pStr = '_wLoHi';
            elseif pCt == 2
                pStr = '_NoLoHi';
            end

            [~,~,VesselDiameterHistBin] = histcounts(cat(1,GroupID_AllVidData_gCell_In{~FlagID(2:end),3}),UI_In.VesselDiameterHistEdge_um);
            cmap = colormap(jet(numel(UI_In.VesselDiameterHistEdge_um)-1)); % Rainbow; Small diameter = Blue; Large diameter = red
            cmap = vertcat(cmap,[0.8 0.8 0.8]); 
            VesselDiameterHistBin(VesselDiameterHistBin==0) = numel(UI_In.VesselDiameterHistEdge_um); 
            cmap = cmap(VesselDiameterHistBin,:);
            close all;

            FigHandle = figure('Position',[10 50 1500 800],'visible','off');
            tlo = tiledlayout(1,1,'TileSpacing','Compact','Padding','tight');

            ax = nexttile(tlo);
            hold(ax, 'on');
            grid on;

            plot(GroupID_CumTorVvsCumConc_Trim(:,1,1),GroupID_CumTorVvsCumConc_Trim(:,2:end,1),'LineWidth',1.5); % Mean Cum WBC Conc
            if pCt == 1
                plot(GroupID_CumTorVvsCumConc_Trim(:,1,2),GroupID_CumTorVvsCumConc_Trim(:,2:end,2),'--','LineWidth',0.5); % Lo Cum WBC Conc
                plot(GroupID_CumTorVvsCumConc_Trim(:,1,3),GroupID_CumTorVvsCumConc_Trim(:,2:end,3),'--','LineWidth',0.5); % Hi Cum WBC Conc
            end
            colororder(cmap);
            legend(GroupID_AllVidData_gCell_In{~FlagID(2:end),1,1},'Orientation','vertical','Location','eastoutside');
            xlabel(strcat('Cum',xStr),'Interpreter','none');
            ylabel('Cum WBC Conc (K/uL)','Interpreter','none');
            xlim([0 CumTorVLim]);
            xticks(0:10:CumTorVLim);
            ylim([0.5*floor(min(GroupID_CumTorVvsCumConc_Trim(:,2:end,2),[],'all')/0.5) 0.5*ceil(max(cellfun(@(x) max(x(:,14),[],1), GroupID_AllVidData_gCell_In(:,2)),[],1)/0.5)]);
            yticks([0.5*floor(min(GroupID_CumTorVvsCumConc_Trim(:,2:end,2),[],'all')/0.5):0.5:0.5*ceil(max(cellfun(@(x) max(x(:,14),[],1), GroupID_AllVidData_gCell_In(:,2)),[],1)/0.5)]);
            set(gca,'TickDir','out');

            % = Save Plot
            Plot_GroupID_CumTorVvsCumConc_lpmCell_Out{torvlCt,pCt,mCt} = export_fig('-r150');
            PlotLabel_GroupID_CumTorVvsCumConc_lpmCell_Out{torvlCt,pCt,mCt} = strcat(NormStr,MethodStr,pStr,'_0-',num2str(CumTorVLim));

            close all;
            clearvars FigHandle ax cmap tlo;

        end
        
        % === Organize plot data for Saving
        WinDiameterCBCTbl = table;
        WinDiameterCBCTbl.GroupID = GroupID_AllVidData_gCell_In(~FlagID(2:end),1);
        WinDiameterCBCTbl.WinVesselDiameter_um = cell2mat(GroupID_AllVidData_gCell_In(~FlagID(2:end),3));
        WinDiameterCBCTbl.WBCConcREF_KuL = cell2mat(GroupID_AllVidData_gCell_In(~FlagID(2:end),14));
        if NormOption==0 | NormOption==1 | NormOption==2 
            WinDiameterCBCTbl.WBCConcCumTorVLim_Mean_KuL = GroupID_CumTorVvsCumConc_Trim(end,2:end,1)';
            WinDiameterCBCTbl.WBCConcCumTorVLim_Lo_KuL = GroupID_CumTorVvsCumConc_Trim(end,2:end,2)';
            WinDiameterCBCTbl.WBCConcCumTorVLim_Hi_KuL = GroupID_CumTorVvsCumConc_Trim(end,2:end,3)';
        elseif NormOption==3 | NormOption==4 
            WinDiameterCBCTbl.WBCConcCumTorVLim_Mean_KuL = (arrayfun(@(x) GroupID_CumTorVvsCumConc_Trim(CumConc_LastNonNaNIdx(x),x,1),2:width(GroupID_CumTorVvsCumConc_Trim)))';
            WinDiameterCBCTbl.WBCConcCumTorVLim_Lo_KuL = (arrayfun(@(x) GroupID_CumTorVvsCumConc_Trim(CumConc_LastNonNaNIdx(x),x,2),2:width(GroupID_CumTorVvsCumConc_Trim)))';
            WinDiameterCBCTbl.WBCConcCumTorVLim_Hi_KuL = (arrayfun(@(x) GroupID_CumTorVvsCumConc_Trim(CumConc_LastNonNaNIdx(x),x,3),2:width(GroupID_CumTorVvsCumConc_Trim)))';
        end

         % = Prepare "GroupID_CumTorVvsCumConc_Trim_Tbl " for
         % Excel worksheets saving 
        GroupID_Label = GroupID_AllVidData_gCell_In(~FlagID(2:end),1);
        GroupID_CumTorVvsCumConc_Trim_Mean_Tbl = ...
            array2table(GroupID_CumTorVvsCumConc_Trim(:,:,1),'VariableNames',vertcat(strcat('Cum',xStr),cellfun(@(x) strcat(x,MethodStr), GroupID_Label,'UniformOutput',false)));
        GroupID_CumTorVvsCumConc_Trim_Lo_Tbl = ...
            array2table(GroupID_CumTorVvsCumConc_Trim(:,2:end,2),'VariableNames',vertcat(cellfun(@(x) strcat(x,MethodStr,'_Lo'), GroupID_Label,'UniformOutput',false)));
        GroupID_CumTorVvsCumConc_Trim_Hi_Tbl = ...
            array2table(GroupID_CumTorVvsCumConc_Trim(:,2:end,3),'VariableNames',vertcat(cellfun(@(x) strcat(x,MethodStr,'_Hi'), GroupID_Label,'UniformOutput',false)));
        
        Tbl_GroupID_CumTorVvsCumConc_lmCell_Out{torvlCt,1,mCt} = ...
            horzcat(GroupID_CumTorVvsCumConc_Trim_Mean_Tbl,GroupID_CumTorVvsCumConc_Trim_Lo_Tbl,GroupID_CumTorVvsCumConc_Trim_Hi_Tbl);
        TblLabel_GroupID_CumTorVvsCumConc_lmCell_Out{torvlCt,1,mCt} = strcat(NormStr,MethodStr,'_0-',num2str(CumTorVLim));
        
        clearvars GroupID_CumTorVvsCumConc_Trim_Mean_Tbl GroupID_CumTorVvsCumConc_Trim_Lo_Tbl GroupID_CumTorVvsCumConc_Trim_Hi_Tbl; 

        % = Load PlotData into output cell
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,1} = CumTorVLim;
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,2} = NormOption;
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,3} = GroupID_CumTorVvsCumConc_Trim; % 20240120: Matrix of depth=3
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,4} = WinDiameterCBCTbl;
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,5} = UI_In.ASCOption;
        if UI_In.ASCOption == 1 
            GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,6} = UI_In.ASC.IncrVolBin.ASCRatioDiffLim;
        elseif UI_In.ASCOption == 2 
            GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,6} = UI_In.ASC.CumVolBin.ASCRatioDiffLim;
        end
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,7} = GroupID_CumTorVvsCumConc_Trim; 
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,8} = WinDiameterCBCTbl; 
        
        % clearvars CumTorVLim FlagID VesselDiameterHistBin cmap;
        % clearvars GroupID_CumTorVvsCumConc_Trim WinDiameterCBCTbl;
        % clearvars GroupID_Label GroupID_CumTorVvsCumConc_Trim_Tbl GroupID_CumTorVvsCumConc_Trim_Mean_Tbl GroupID_CumTorVvsCumConc_Trim_Lo_Tbl GroupID_CumTorVvsCumConc_Trim_Hi_Tbl;

    end

    GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_Out{mCt,1} = GroupID_CumTorVvsCumConc_PlotData_torvlCell;

end


close all;
clearvars FigHandle ax NormOption NormStr tlCt CumTorVLim_All ugCt2 pCt mCt MethodStr xStr;
clearvars GroupID_CumTorVvsCumConc GroupID_CumTorVvsCumConc_PlotData_torvlCell;

