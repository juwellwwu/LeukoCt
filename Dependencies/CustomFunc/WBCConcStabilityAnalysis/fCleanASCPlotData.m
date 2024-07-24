function [GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_Out,...
    Tbl_GroupID_CumTorVvsCumConc_Stable_lmCell_Out, TblLabel_GroupID_CumTorVvsCumConc_Stable_lmCell_Out] = ...
    fCleanASCPlotData(GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_In,GroupID_AllVidData_WBCConcStability_Tbl_In,UI_In,TimeorVolOption_In)

fprintf('Keep only ASC-Stable, ImgTime (or Flow Volume) > CumTovVLim Group IDs in GroupID_CumTorVvsCumConc_PlotData_torvlCell ...\n');


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


%% Clean up

GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_Out = GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_In;

% = Cell for holding Excel CumTorV vs CumConc Table
Tbl_GroupID_CumTorVvsCumConc_Stable_lmCell_Out = cell(numel(CumTorVLim_All),1,2); 
TblLabel_GroupID_CumTorVvsCumConc_Stable_lmCell_Out = cell(numel(CumTorVLim_All),1,2); 

for mCt = 1:2 

    if mCt == 1
        MethodStr = '_OB';
    elseif mCt == 2
        MethodStr = '_SK';
    end

    GroupID_CumTorVvsCumConc_PlotData_torvlCell = GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_In{mCt,1};

    for torvlCt = 1:height(GroupID_CumTorVvsCumConc_PlotData_torvlCell) 

        CumTorVLim = GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,1};

        % = Determine GroupIDs with ImgTime (s) -OR- FlowVol (nL) > CumTorVLim & pass ASC stability
        % criteria
        if TimeorVolOption_In == 1 
            GroupID_LongStable = ...
                intersect(GroupID_CumTorVvsCumConc_PlotData_torvlCell{1,8}.GroupID,... 
                GroupID_AllVidData_WBCConcStability_Tbl_In.GroupID(GroupID_AllVidData_WBCConcStability_Tbl_In.ASCImgTimeStart_s>eps)); 
        elseif TimeorVolOption_In == 2 
            GroupID_LongStable = ...
                intersect(GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,8}.GroupID,... 
                GroupID_AllVidData_WBCConcStability_Tbl_In.GroupID(GroupID_AllVidData_WBCConcStability_Tbl_In.ASCCumFlowVolStart_nL>eps)); 
        end
        FlagID = ~ismember(GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,8}.GroupID,GroupID_LongStable);

        % = Remove GroupIDs in Col 7 (GroupID_CumTorVvsCumConc_Trim), Col 8 (WinDiameterCBCTbl) that do not meet ASC WBC Conc Stability criteria
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,7}(:,find(FlagID)+1,:) = []; 
        GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,8}(FlagID,:) = [];

        % = Prepare Col 7-8 data tables for Saving
        GroupID_Label = GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,8}.GroupID;

        GroupID_CumTorVvsCumConc_Trim_Stable_Mean_Tbl = ...
            array2table(GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,7}(:,:,1),'VariableNames',vertcat('CumImgTime_s',cellfun(@(x) strcat(x,'_Mean'), GroupID_Label,'UniformOutput',false)));
        GroupID_CumTorVvsCumConc_Trim_Stable_Lo_Tbl = ...
            array2table(GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,7}(:,2:end,2),'VariableNames',vertcat(cellfun(@(x) strcat(x,'_Lo'), GroupID_Label,'UniformOutput',false)));
        GroupID_CumTorVvsCumConc_Trim_Stable_Hi_Tbl = ...
            array2table(GroupID_CumTorVvsCumConc_PlotData_torvlCell{torvlCt,7}(:,2:end,3),'VariableNames',vertcat(cellfun(@(x) strcat(x,'_Hi'), GroupID_Label,'UniformOutput',false)));
        
        Tbl_GroupID_CumTorVvsCumConc_Stable_lmCell_Out{torvlCt,1,mCt} = ...
            horzcat(GroupID_CumTorVvsCumConc_Trim_Stable_Mean_Tbl,GroupID_CumTorVvsCumConc_Trim_Stable_Lo_Tbl,GroupID_CumTorVvsCumConc_Trim_Stable_Hi_Tbl);
        TblLabel_GroupID_CumTorVvsCumConc_Stable_lmCell_Out{torvlCt,1,mCt} = strcat(NormStr,MethodStr,'_0-',num2str(CumTorVLim));
        
        clearvars GroupID_CumTorVvsCumConc_Trim_Stable_Mean_Tbl GroupID_CumTorVvsCumConc_Trim_Stable_Lo_Tbl GroupID_CumTorVvsCumConc_Trim_Stable_Hi_Tbl; 
        clearvars GroupID_LongStable FlagID;
 
    end

    GroupID_CumTorVvsCumConc_PlotData_torvlCell_mCell_Out{mCt,1} = GroupID_CumTorVvsCumConc_PlotData_torvlCell;

    clearvars GroupID_CumTorVvsCumConc_PlotData_torvlCell;

end

clearvars NormOption NormStr mCt MethodStr torvlCt;

