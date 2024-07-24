clc; close all hidden; fclose('all'); clearvars;


%% USER INPUT

% = Cell Counting Pipeline Output Folder *****
UI.PipeOut_FoldernameString = '/Users/j/Documents/MATLAB/LeukoCt/Data_In/Stability_In/'; % End w/ "/"

% = Output Root Folder (include / at end) *****
UI.Output_RootFoldernameString = '/Users/j/Documents/MATLAB/LeukoCt/Data_Out/Stability/';

% = WBC Stable Concentration Criteria
UI.ASCOption = 2; % Do not change

UI.ASC.CumVolBin.ASCRatioDiffLim = 0.05; % default: 0.05
UI.ASC.CumVolBin.ReportIncVol_Skip = 1E-6;  % default: 1E-6 
UI.ASC.CumVolBin.ReportIncVol_Step = [2.0:0.01:4.0]; % default: [2.0:0.01:4.0]
UI.ASC.CumVolBin.MinStableBinCt = 5; % default: 5
UI.ASC.CumVolBin.MinStableVol = 20; % default: 20
UI.ASC.CumVolBin.MinStableTime = -9999; % default: -9999 (off)

% % = Interpotation Step Size
UI.CumVolIntp_Step = 0.01; % default: 0.01 
UI.CumTimeIntp_Step = 0.01; % default: 0.01

% 

% = How flagged tMidPt are treated 
% FlagtMidPtCutOption:
% 0 = Do not remove Flagged tMidPt of any flags
% 1 = Remove tMidPt flagged by VctyEr_OB only
UI.FlagtMidPtCutOption = 1; % Do not change

% = Option to skip individual GroupID plots
UI.GroupIDPlotOption = 'y';


%% (The following user-input does not affect reproducibility of published data) 

% = How CellCt Window's Vessel Diameter (um) is grouped for statistics 
UI.VesselDiameterHistEdge_um = [10.0 16.0 22.0 28.0]';

% = WBC Count Correction
% 0 to use original auto WBC count, 1 to subtract expected false + count
% based on # video frames and threshold ScMAD multiple ("multScMAD_CellCt")
UI.CtCorrOption = 0; % default: 0

% = Cumulative ImgTime or FlowVol Normalization Option
% 0: No Normalization, GroupID w/ CumTime or CumVol >CumTimeLim or CumVolLim
% 1: Normalization Method 1: Divide by largest value designated by
% CumTimeLim_ALL or CumVolLim_ALL; GroupID w/ CumTime or CumVol > CumTimeLim or CumVolLim 
% 2: Normalization Method 2: Divide by largest value designated by
% CumTimeLim_ALL or CumVolLim_ALL to get x, then abs(x-1). All resulting values are positive; 
% GroupID w/ CumTime or CumVol >CumTimeLim or CumVolLim  
% 3: Same as Normalization Method 0 (No Normalization), All GroupID
% 4: Same as Normalization Method 1, All GroupID
UI.CumTimeVolvsCumConcNormOption = 3; 
UI.CumTimeLim_All = 900; % One or more Max ImgTime (s) to include in CumTime vs CumWBCConc Plot; ex. [60:60:180] 
UI.CumVolLim_All = 550; % One or more Max Cum Flow Vol (nL) to include in CumVol vs CumWBCConc Plot; ex. [50:50:150] 


%% Dummy input variables [Ignore values]

UI.ASC.IncrVolBin.ASCRatioDiffLim = UI.ASC.CumVolBin.ASCRatioDiffLim;
UI.ASC.IncrVolBin.ReportIncVol_Skip = UI.ASC.CumVolBin.ReportIncVol_Skip;
UI.ASC.IncrVolBin.ReportIncVol_Step = UI.ASC.CumVolBin.ReportIncVol_Step;
UI.ASC.IncrVolBin.MinStableBinCt = UI.ASC.CumVolBin.MinStableBinCt;
UI.ASC.IncrVolBin.MinStableVol = UI.ASC.CumVolBin.MinStableVol;
UI.ASC.IncrVolBin.MinStableTime = UI.ASC.CumVolBin.MinStableTime;

UI.FlagtMidPtCut.CellCtDiameterLim_um = 9999;

UI.CtCorrOption = 0;


%% Prepare Output Folder

OutputFolder = strcat(UI.Output_RootFoldernameString,'WBCConcStabilityAnalysis'); % No '/' at end

timestamp = string(datetime('now'),'yyMMddHHmm');
SaveFilePath = strcat(OutputFolder,'_',timestamp,'/');
mkdir(SaveFilePath);

% % % SaveTBinPlotFilePath = strcat(SaveFilePath,'Plot_TimeBin/'); % One Plot per GroupID
% % % if exist(SaveTBinPlotFilePath,'dir')==7
% % %     rmdir(SaveTBinPlotFilePath,'s'); 
% % % end
% % % mkdir(SaveTBinPlotFilePath);

% % % SaveTBinAllPlotFilePath = strcat(SaveFilePath,'Plot_TimeBin_AllGroupID/');
% % % if exist(SaveTBinAllPlotFilePath,'dir')==7
% % %     rmdir(SaveTBinAllPlotFilePath,'s');
% % % end
% % % mkdir(SaveTBinAllPlotFilePath);

if UI.ASCOption==2 % Apply Stability criteria to Cumulative Flow Vol Bin
    SaveVBinPlotFilePath = strcat(SaveFilePath,'Plot_FlowVolBin_CumWBCConc/');
end
if exist(SaveVBinPlotFilePath,'dir')==7
    rmdir(SaveVBinPlotFilePath,'s'); 
end
mkdir(SaveVBinPlotFilePath);


%% Identify GroupID for each directory
% Each GroupID present 1 vessel, multiple vessels may exist in 1 video

fprintf('Identify GroupID (Vessel identifier) for all videos in UI.PipeOut_FoldernameString directory ...\n');

% = Determine list of AUTO_PIPELINE directories in User Input
PipeOut_FolderFileList = dir(fullfile(UI.PipeOut_FoldernameString));
PipeOut_FolderFileList(~cat(1,PipeOut_FolderFileList.isdir)) = []; % Keep only directories
PipeOut_FolderFileList(ismember({PipeOut_FolderFileList.name},{'.','..'})) = [];
PipeOut_FolderFileList(find(cellfun(@(x) isempty(x), regexp({PipeOut_FolderFileList.name},'\d+_P\d+_V\d+')))) = []; % Remove non-sample directories

% = Prepare Cell summarizing each video (multiple per GroupID)
% "GroupID_Data_fCell" will hold GroupID (Vessel Identifier), VideoID in each directory
% 1) Group ID (Vessel Identifier)
% 2) Video ID
% 3) VideoLength (s) of video
% 4) Blood Flow Volume (nL) of video, OB (variation from orient band slope angles)
% 5) Blood Flow Volume (nL) of video, SK (variation from skeletons)
% 6) Window diameter (um) of video

GroupID_Data_fCell = repmat({'','',[],[],[],[]},height(PipeOut_FolderFileList),1); 

for fCt = 1:height(PipeOut_FolderFileList)

    SampleIDStr_Temp = PipeOut_FolderFileList(fCt).name;

    % = Load Group ID = Patient ID + VesselID as class 'character' 
    expr_PatientID = 'P\d+';
    [StartIdx_PatientID,EndIdx_PatientID] = regexp(SampleIDStr_Temp,expr_PatientID);
    if ~isempty(StartIdx_PatientID)
        PatientIDString = convertCharsToStrings(SampleIDStr_Temp(StartIdx_PatientID+1:EndIdx_PatientID));
    else
        PatientIDString = "";
    end

    expr_VesselID = 'V\d+';
    expr_AppendSampleIDStr = '_[a-z]'; % 20240105: specify which vessel within one video
    [StartIdx_VesselID,EndIdx_VesselID] = regexp(SampleIDStr_Temp,expr_VesselID);
    if ~isempty(StartIdx_VesselID)
        VesselIDString = convertCharsToStrings(SampleIDStr_Temp(StartIdx_VesselID+1:EndIdx_VesselID)); % convert to string
        [StartIdx_AppendSampleIDStr,EndIdx_AppendSampleIDStr] = regexp(SampleIDStr_Temp,expr_AppendSampleIDStr);
        if ~isempty(StartIdx_AppendSampleIDStr)
            VesselIDString = strcat(VesselIDString,convertCharsToStrings(SampleIDStr_Temp(StartIdx_AppendSampleIDStr+1:EndIdx_AppendSampleIDStr)));
        end
    else
        VesselIDString = "";
    end

    expr_VesselVideoID = 'V\d+_\d+';
    expr_VideoID = '_\d+'; 
    [StartIdx_VesselVideoID,EndIdx_VesselVideoID] = regexp(SampleIDStr_Temp,expr_VesselVideoID);
    if ~isempty(StartIdx_VesselVideoID)
        VesselVideoIDString = SampleIDStr_Temp(StartIdx_VesselVideoID:EndIdx_VesselVideoID); % convert to string
        [StartIdx_VideoIDStr,EndIdx_VideoIDStr] = regexp(VesselVideoIDString,expr_VideoID);
        if ~isempty(StartIdx_VideoIDStr)
            VideoIDString = convertCharsToStrings(VesselVideoIDString(StartIdx_VideoIDStr+1:EndIdx_VideoIDStr));
        end
    else
        VideoIDString = "";
    end

    GroupID_Data_fCell{fCt,1} = convertStringsToChars(strcat(PatientIDString,VesselIDString));
    GroupID_Data_fCell{fCt,2} = convertStringsToChars(VideoIDString);

    clearvars StartIdx* EndIdx* expr_* SampleIDStr_Temp PatientIDString VesselIDString VideoIDString VesselVideoIDString;

end

clearvars fCt;


%% Prepare Time-Binned, cumulative flow vol and WBC conc data for each GroupID (Vessel identifier) 

fprintf('Prepare time-binned, cumulative flow vol and WBC conc data for each GroupID (Vessel identifier) ...\n');

% = Find unique GroupIDs
% Will group all videos belong to each Group ID for cumulative blood vol
% and WBC conc analysis
[GrpID_Idx,GrpID_uq] = findgroups(GroupID_Data_fCell(:,1));

% 'GroupID_AllVidData_gCell' has following columns:
% Col 1) GroupID
% Col 2) tMidPt Data (Volume, Conc, Flag etc) from PIPELINE
% Col 3) Mean Window Vessel Diameter (um) for GroupID
% Col 4) WBC Conc & Imaging Time Data based on equal Incremental Blood Flow Vol 
% (RData)
% Col 5) Auto Stable Criteria (ASCParam)
% Col 6) WBC location (sec; ref to start video 1) + WBC Z Score 
% Col 7) Stable WBC Conc: Imaging Time (s) Start
% Col 8) Stable WBC Conc: Cum Flow Vol (nL) Start
% Col 9) Stable WBC Conc: Best Flow Vol Bin Size (nL)
% Col 10) Stable WBC Conc: Start Flow Vol Bin Idx
% Col 11) Stable WBC Conc: End Flow Vol Bin Idx
% Col 12) Stable WBC Conc: Max WBC Conc Ratio-Difference-from-1 w/in Stable Concentration volume 
% Col 13) Stable WBC Conc: Final (Incremental or Cumulative) Flow Vol Bin
% WBC Concentration (K/uL)
% Col 14) REF WBC Concentration (K/uL) from CBC

GroupID_AllVidData_gCell = repmat({'',{},zeros(1,1,'single'),{},{},{}, ... % Col 1-6
    -9999*ones(1,1,'single'),-9999*ones(1,1,'single'),-9999*ones(1,1,'single'),-9999*ones(1,1,'single'),-9999*ones(1,1,'single'),-9999*ones(1,1,'single'),... % Col 7-12
    -9999*ones(1,1,'single'),-9999*ones(1,1,'single')},... % Col 13-
    numel(GrpID_uq),1); % Size of repmat

for ugCt = 1:numel(GrpID_uq) 

    fprintf(strcat('GroupID:',32,num2str(ugCt),'/',num2str(numel(GrpID_uq)),'...\n'));

    GroupID_AllVidData_gCell{ugCt,1} = GrpID_uq{ugCt,1}; % Col 1: Group ID   
   
    vList_Temp = PipeOut_FolderFileList(GrpID_Idx==ugCt); % folder list for unique GroupID; each folder = each video

    % = Prepare data from all videos of the vessel
  
    [tMidPtData,CellLocZScore_Cell,RefConc_KuL] = fCombineVesselVidData(vList_Temp,4,UI);
    
    % = Load per video data into "GroupID_Data_fCell":
    [vidGrp,vidGrpID] = findgroups(tMidPtData(:,1)); 
    vidRowIdx = find(GrpID_Idx==ugCt); 
    for vCt2 = 1:numel(vidRowIdx)
        GroupID_Data_fCell{vidRowIdx(vCt2),3} = max(tMidPtData(tMidPtData(:,1)==vCt2,3),[],1); 
        GroupID_Data_fCell{vidRowIdx(vCt2),4} = max(tMidPtData(tMidPtData(:,1)==vCt2,31),[],1); 
        GroupID_Data_fCell{vidRowIdx(vCt2),5} = max(tMidPtData(tMidPtData(:,1)==vCt2,34),[],1); 
        GroupID_Data_fCell{vidRowIdx(vCt2),6} = max(tMidPtData(tMidPtData(:,1)==vCt2,18),[],1); 
    end
    clearvars vCt2 vidGrp vidGrpID;

    % = Load all-video data in 'GroupID_AllVidData_gCell':
    GroupID_AllVidData_gCell{ugCt,2} = tMidPtData;
    GroupID_AllVidData_gCell{ugCt,3} = ...
         sum(cat(1,GroupID_Data_fCell{vidRowIdx,6}).*cat(1,GroupID_Data_fCell{vidRowIdx,4}),1)./sum(cat(1,GroupID_Data_fCell{vidRowIdx,4}),1);
    GroupID_AllVidData_gCell{ugCt,6} = CellLocZScore_Cell;
    GroupID_AllVidData_gCell{ugCt,14} = RefConc_KuL;
    clearvars tMidPtData CellLocZScore_Cell RefConc_KuL vidRowIdx;
   
end

clearvars GrpID_Idx GrpID_uq vList_Temp PipeOut_FolderFileList;


%% Plot tMidPt Data (all videos combined for each vessel)

% % % if UI.GroupIDPlotOption == 'y'
% % % 
% % %     [Fig_VidCombotMidPtData_Cell] = fPlotVidCombotMidPtData(GroupID_AllVidData_gCell);
% % % 
% % %     for ugCt = 1:height(GroupID_AllVidData_gCell)
% % %         if ~isempty(Fig_VidCombotMidPtData_Cell{ugCt,1})
% % %             imwrite(Fig_VidCombotMidPtData_Cell{ugCt,1},strcat(SaveTBinPlotFilePath,'Plot_TimeBin_CumFlowVolImgTimeWBCConc_', GroupID_AllVidData_gCell{ugCt,1},'_',timestamp,'.tif'));
% % %         end
% % %     end
% % % 
% % % end

clearvars ugCt Fig_VidCombotMidPtData_Cell;


%% Interpolate Cumulative WBC Concentration (K/uL) vs imaging time (s) for each GroupID (Vessel identifier)

[GroupID_CumTimevsCumConc_mCell] = fInterpolateCumWBCConc(GroupID_AllVidData_gCell,1,UI);


%% Interpolate Cum WBC Concentration (K/uL) vs Cum Blood Vol (nL) for each GroupID (Vessel Identifier)

[GroupID_CumVolvsCumConc_mCell] = fInterpolateCumWBCConc(GroupID_AllVidData_gCell,2,UI);


%% Plot: Cum Imaging Time (sec) (x) vs WBC Concentration (K/uL) (y) for each GroupID (Vessel Identifier)
% Cum WBC Conc (y) vs Cum Img Time (x) line plot of all cum img time > CumTimeLim Group IDs.
% May be Stable or Unstable

% = Prepare Plot Save Directory
% % % SaveTBinAllPlotCumTimevsCumConcFilePath = ...
% % %     strcat(SaveTBinAllPlotFilePath,'CumImgTime_vs_CumWBCConc',strcat('_Norm',num2str(UI.CumTimeVolvsCumConcNormOption)),'/'); 
% % % if exist(SaveTBinAllPlotCumTimevsCumConcFilePath,'dir')==7
% % %     rmdir(SaveTBinAllPlotCumTimevsCumConcFilePath,'s'); % important
% % % end
% % % mkdir(SaveTBinAllPlotCumTimevsCumConcFilePath);

% = Plot and Prepare Plot Data
% "GroupID_CumTorVvsCumConc_PlotData_torvlCell " cell holds plot data
% "TorV": T = Cumulative Imaging Time (s); V = Cumulative Flow Volume (nL)
% # Rows = # CumTorVLim
% Col 1) CumTorVLim
% Col 2) Normalization Method (0-4)
% Col 3) GroupID_CumTorVvsCumConc Data for GroupID w/ ImgTime -OR- Flow Vol >=
% CumTorVLimit (= GroupID_CumTorVvsCumConc_Trim)
% Col 4) Window Vessel Diameter (um), REFCellConc (K/uL), for GroupID w/ ImgTime -OR- Flow Vol >=
% CumTorVLimit (= WinDiameterCBCTbl)
% Col 5) ASC WBC Conc Stability Option: 2 if use Cumulative Vol Bin (= UI.ASCOption)
% Col 6) ASC WBC Conc Stability: Max allowed WBC Conc
% Ratio-Difference-from-1 (= UI.ASCRatioDiffLim_CumVolBin)
% Col 7) GroupID_CumTorVvsCumConc Data for GroupID w/ ImgTime -OR- Flow Vol >=
% CumTorVLimit & Pass ASC Stabillty
% Col 8) Window Vessel Diameter (um), REFCellConc (K/uL), for GroupID w/ ImgTime -OR- Flow Vol >=
% CumTorVLimit (= WinDiameterCBCTbl) & Pass ASC Stabillty

[GroupID_CumTimevsCumConc_PlotData_tlCell_mCell, ...
    Plot_GroupID_CumTimevsCumConc_lpmCell, PlotLabel_GroupID_CumTimevsCumConc_lpmCell, ...
    Tbl_GroupID_CumTimevsCumConc_lmCell, TblLabel_GroupID_CumTimevsCumConc_lmCell] = ...
    fPlotCumWBCConc(GroupID_CumTimevsCumConc_mCell,GroupID_AllVidData_gCell,UI,1); % Vector input = How to divide Window Vessel Diameter into groups

% Save plot
% % % for torvlCt = 1:height(Plot_GroupID_CumTimevsCumConc_lpmCell) % # Cum time limit
% % %     for pCt = 1:1 
% % %         for mCt = 1:1 
% % %             imwrite(Plot_GroupID_CumTimevsCumConc_lpmCell{torvlCt,pCt,mCt},strcat(SaveTBinAllPlotCumTimevsCumConcFilePath,'Plot_CumTimevsCumWBCConc',PlotLabel_GroupID_CumTimevsCumConc_lpmCell{torvlCt,pCt,mCt},'s_',timestamp,'.tif'));
% % %         end
% % %     end
% % % end

% Save table
% % % for torvlCt = 1:height(Tbl_GroupID_CumTimevsCumConc_lmCell) % # Cum time limit
% % %     for mCt = 1:1 
% % %         writetable(Tbl_GroupID_CumTimevsCumConc_lmCell{torvlCt,1,mCt},strcat(SaveTBinAllPlotCumTimevsCumConcFilePath,'Data_CumTimevsCumWBCConc',TblLabel_GroupID_CumTimevsCumConc_lmCell{torvlCt,1,mCt},'s_',timestamp,'.csv'))
% % %     end
% % % end

clearvars Plot*_GroupID_CumTimevsCumConc_lpmCell Tbl*_GroupID_CumTimevsCumConc_lmCell torvlCt pCt mCt;


%% Plot: Cum Blood Vol (nL) (x) vs Cum WBC Concentration (K/uL) (y) for each GroupID (Vessel Identifier)
% Cum WBC Conc (y) vs Cum Flow Vol (x) line plot of all cum flow vol > CumVolLim Group IDs.
% May be ASC-Stable or Unstable.

% = Prepare Plot Save Directory
% % % SaveTBinAllPlotCumVolvsCumConcFilePath = strcat(SaveTBinAllPlotFilePath,'CumFlowVol_vs_CumWBCConc',strcat('_Norm',num2str(UI.CumTimeVolvsCumConcNormOption)),'/'); 
% % % if exist(SaveTBinAllPlotCumVolvsCumConcFilePath,'dir')==7
% % %     rmdir(SaveTBinAllPlotCumVolvsCumConcFilePath,'s'); % important
% % % end
% % % mkdir(SaveTBinAllPlotCumVolvsCumConcFilePath);

% = Plot
[GroupID_CumVolvsCumConc_PlotData_vlCell_mCell, ...
    Plot_GroupID_CumVolvsCumConc_lpmCell, PlotLabel_GroupID_CumVolvsCumConc_lpmCell, ...
    Tbl_GroupID_CumVolvsCumConc_lmCell, TblLabel_GroupID_CumVolvsCumConc_lmCell] = ...
    fPlotCumWBCConc(GroupID_CumVolvsCumConc_mCell,GroupID_AllVidData_gCell,UI,2); % Vector input = How to divide Window Vessel Diameter into groups

% Save plot
% % % for torvlCt = 1:height(Plot_GroupID_CumVolvsCumConc_lpmCell) % # Cum volume limit
% % %     for pCt = 1:1 
% % %         for mCt = 1:1 
% % %             imwrite(Plot_GroupID_CumVolvsCumConc_lpmCell{torvlCt,pCt,mCt},strcat(SaveTBinAllPlotCumVolvsCumConcFilePath,'Plot_CumVolvsCumWBCConc',PlotLabel_GroupID_CumVolvsCumConc_lpmCell{torvlCt,pCt,mCt},'nL_',timestamp,'.tif'));
% % %         end
% % %     end
% % % end

% Save table
% % % for torvlCt = 1:height(Tbl_GroupID_CumVolvsCumConc_lmCell) % # Cum volume limit
% % %     for mCt = 1:1 
% % %         writetable(Tbl_GroupID_CumVolvsCumConc_lmCell{torvlCt,1,mCt},strcat(SaveTBinAllPlotCumVolvsCumConcFilePath,'Data_CumVolvsCumWBCConc',TblLabel_GroupID_CumVolvsCumConc_lmCell{torvlCt,1,mCt},'nL_',timestamp,'.csv'))
% % %     end
% % % end

clearvars Plot*_GroupID_CumVolvsCumConc_lpmCell Tbl*_GroupID_CumVolvsCumConc_lmCell torvlCt pCt mCt;


%% ===== Start ASC (Auto Stability Criteria) Analysis

%% Interpolate Cumulative Blood Flow Vol (nL) (x; Col1) vs Imaging time (s) (y;Col2-end) for each GroupID (Vessel identifier)

[GroupID_CumVolvsCumTime_mCell] = fInterpolateCumImgTime(GroupID_AllVidData_gCell,UI);


%% Calculate Starting Time and Cumulative Flow Vol at which WBC Conc achieves Stability

fprintf('Calculating Starting Time and Cumulative Flow Vol at which WBC Conc achieves Stability ...\n');

GroupID_CumVolvsCumTime = GroupID_CumVolvsCumTime_mCell{1,1}; % Method OB data

% = Test Parameters
if UI.ASCOption == 2 % Apply Stability Criteria to Cumulative Flow Vol Bin
    ASC_ReportIncVol_Skip = UI.ASC.CumVolBin.ReportIncVol_Skip;
    ASC_ReportIncVol_Step = UI.ASC.CumVolBin.ReportIncVol_Step;
    ASCRatioDiffLim = UI.ASC.CumVolBin.ASCRatioDiffLim;
end

for ugCt = 1:height(GroupID_AllVidData_gCell)

    % ASCParam = GroupID_AllVidData_gCell{ugCt,5} column information
    % Col 1). ASC_ReportIncVol_Step
    % Col 2). # Non-Skip, Complete (Full Volume) Incremental Flow Volume Bin
    % Col 3). # Complete Inc Flow Vol Bins w/ stable WBC Conc (-9999 if
    % none; ie, total recorded Flow Vol < ASC_ReportIncVol_Skip+ASC_ReportIncVol_Step)
    % Col 4). Start of Bin Idx of Complete Inc Flow Vol Bins w/ stable WBC Conc 
    % Col 5). End of Bin Idx of Complete Inc Flow Vol Bins w/ stable WBC Conc 
    % Col 6). Start of CumFlowVol (nL) w/ stable WBC Conc 
    % Col 7). End of CumFlowVol (nL) w/ stable WBC Conc
    % Col 8). FlowVol (nL) w/ stable WBC Conc (= Col5-Col4)
    % Col 9). Max Ratio-Difference-from-1 found while fullfilling stable criteria
    % Col 10). BestVolStep_ASCParamRow: 1 if row is best row, 0 otherwise

    ASCParam = -9999*ones(numel(ASC_ReportIncVol_Step),10,'single');
    ASCParam(:,9) = 9999;

    GroupID_AllVidData_gCell{ugCt,4} = cell(1,1,numel(ASC_ReportIncVol_Step)); % For storing RData; 1 RData for each ASC_ReportIncVol_Step

    for stCt = 1:numel(ASC_ReportIncVol_Step) 

        ASCParam(stCt,1) = ASC_ReportIncVol_Step(stCt);

        % = Variable Skip Flow Vol Bin (1st Bin) Size
        CumVolHistEdge = ...
            fliplr(horzcat(GroupID_CumVolvsCumTime(find(~isnan(GroupID_CumVolvsCumTime(:,ugCt+1)),1,'last'),1):-ASC_ReportIncVol_Step(stCt):ASC_ReportIncVol_Skip,0));
        CumVolHistEdge = CumVolHistEdge + UI.CumVolIntp_Step;
        CumVolHistEdge(1) = 0;

        if all(CumVolHistEdge<eps)  
            continue
        end

        [~,~,ReportIncVolHistBin] = ... 
            histcounts(round(GroupID_CumVolvsCumTime(1:find(~isnan(GroupID_CumVolvsCumTime(:,ugCt+1)),1,'last'),1),4), round(CumVolHistEdge,2));

        % = Distribute WBC by report incremental blood volume
        max_wrapper = @(x) max(x,[],'all','omitmissing'); 
        loc_tsMax = splitapply(max_wrapper,GroupID_CumVolvsCumTime(1:find(~isnan(GroupID_CumVolvsCumTime(:,ugCt+1)),1,'last'),ugCt+1),ReportIncVolHistBin);
        loc_tsMax = loc_tsMax(~isnan(loc_tsMax)); 

        loc_tFrHistEdge = vertcat(0,loc_tsMax); 
        
        % = Distribute WBCs into Inc Vol Bins
        [loc_HistCt,~,loc_HistBin] = histcounts(cell2mat(GroupID_AllVidData_gCell{ugCt,6}(:,1)),loc_tFrHistEdge); 

        % = Record # Non-Skip, Complete (Full Volume) Flow Volume Bins
        ASCParam(stCt,2) = numel(loc_tsMax)-1; 

        if ASCParam(stCt,2)<1 
            continue
        end

        % == RData = GroupID_AllVidData_gCell{ugCt,4} 
        % # Planes = # ASC_ReportIncVol_Step (ie, cell within cell)
        % For each RData of each ASC_ReportIncVol_Step:
        % Columns:
        % Col 1) (Report) Incremental FlowVol Bin Idx
        % Col 2) Start of (Report) IncFlowVol (nL, referenced to CumFlowVol) Bin
        % Col 3) End of (Report) IncFlowVol (nL, referenced to CumFlowVol) Bin
        % Col 4) Corresponding Cum Img Time (s) for (Report) IncFlowVol (nL) in Col 1) (Bin Start)
        % Col 5) Corresponding Cum Img Time (s) for (Report) IncFlowVol (nL) in Col 2) (Bin End)
        % Col 6) Incremental WBC Ct in FlowVol between Col 1) and Col 2)
        % Col 7) Cumulative (ASCOption==2) WBC Conc (K/uL) in FlowVol between Col 1) and Col 2)
        % Col 8) (Report) Incremental FlowVol (nL) Bin Center for Plotting
        % Col 9) Stable Criteria: Based on Max Ratio-Difference-from-1 of Cumulative 
        % (ASCOption==2)  WBC Conc to ALL its
        % subsequent incremental or cumulative WBC Concs; if smaller than user-defined number 
        % ("Stable criteria", UI.ASCRatioDiffLim), consider concentration has stabilized
        % Col 10) Bin ICumulative (ASCOption==2) WBC Conc fullfill Stable Criteria? 1 if true, 0 otherwise.
        % Bins that pass may not be continuous, and / or may not include the last 
        % non-skip, complete flow vol bin. 
        % Col 11) Bin Cumulative (ASCOption==2) WBC Conc fullfill Stable Criteria, continuous and 
        % include last non-skip, complete flow vol bin? 1 if true, 0 otherwise.

        RData = zeros(numel(loc_HistCt),11,'single');

        RData(:,1) = (1:1:numel(loc_HistCt))'; 
        RData(:,4) = loc_tFrHistEdge(1:end-1); 
        RData(:,5) = loc_tFrHistEdge(2:end); 
        RData(:,6) = loc_HistCt'; 

        [~,Loc] = ismember(loc_tFrHistEdge(1:end-1),GroupID_CumVolvsCumTime(:,ugCt+1));
        RData(:,2) = GroupID_CumVolvsCumTime(Loc,1); 

        [~,Loc] = ismember(loc_tFrHistEdge(2:end),GroupID_CumVolvsCumTime(:,ugCt+1));
        RData(:,3) = GroupID_CumVolvsCumTime(Loc,1); 

        if UI.ASCOption == 2 
            RData(:,7) = cumsum(RData(:,6))./RData(:,3); 
        end
        RData(:,8) = 0.5*(RData(:,2)+RData(:,3));  

        RData2 = RData;
        RData2(1,:) = []; 

        % = Col 9: Determine Max Ratio-Difference-from-1 
        if UI.ASCOption == 2 
            for cbCt = 1:height(RData2) 
                RatioDiff = abs(RData2(cbCt:end,7)./RData2(end,7)-1);
                RData2(cbCt,9) = max(RatioDiff,[],1);
            end
        end

         % = Col 10: Apply Stable WBC Conc Criteria: Max Ratio-Diff-from-1 < UI.ASCRatioDiffLim 
         RData2(:,10) = RData2(:,9) < ASCRatioDiffLim;

        RData(2:height(RData2)+1,9) = RData2(:,9); 
        RData(2:height(RData2)+1,10) = RData2(:,10); 

        RData(1,9) = -9999; 

        % = Find stretch of 1s including last complete Inc Flow Vol Bin
        rp = regionprops(logical(RData2(:,10)),'Area','PixelIdxList');
        FlagID = zeros(height(rp),1,'logical');
        for rCt= 1:height(rp) 
            FlagID(rCt,1) = ~ismember(height(RData2(:,10)), rp(rCt).PixelIdxList);
        end
        rp(FlagID) = []; 

        if ~isempty(rp)
            ASCParam(stCt,3) = rp.Area; 
            ASCParam(stCt,4) = min(rp.PixelIdxList+1,[],'all'); 
            ASCParam(stCt,5) = max(rp.PixelIdxList+1,[],'all'); 
            if rp.Area>1 
                ASCParam(stCt,6) = RData(rp.PixelIdxList(1)+1,2); 
                ASCParam(stCt,7) = RData(rp.PixelIdxList(end)+1,3); 
                ASCParam(stCt,8) = ASCParam(stCt,7)-ASCParam(stCt,6); 
                ASCParam(stCt,9) = max(RData2(rp.PixelIdxList,9),[],1); 
            end
        end

        ASCParam(stCt,10) = 0; 

        RData(rp.PixelIdxList+1,11) = 1;

        GroupID_AllVidData_gCell{ugCt,4}{1,1,stCt} = RData;

        clearvars ReportIncVolHistBin loc_tsMax loc_HistBin loc_HistCt loc_tFrHistEdge Loc RatioDiff cbCt rp FlagID rCt RData2;

    end

    GroupID_AllVidData_gCell{ugCt,5} = ASCParam;

    clearvars max_wrapper loc_tFrMax Loc CumVolHistEdge;

end

clearvars GroupID_CumVolvsCumTime RData RData2 ASCParam;
clearvars ASC_ReportIncVol_Skip ASC_ReportIncVol_Step ASCRatioDiffLim;
clearvars stCt ugCt;


%% Determine Cum Vol at which WBC concentration is stable for each GroupID

for ugCt = 1:height(GroupID_AllVidData_gCell) 

    if UI.ASCOption == 2 
        BinCtLlim = UI.ASC.CumVolBin.MinStableBinCt;
        MinStableVol = UI.ASC.CumVolBin.MinStableVol;
        MinStableTime = UI.ASC.CumVolBin.MinStableTime;
    end
    
    % Determine which ReportIncVol_Step (=Volume Bin Size) to use
    if UI.ASCOption == 2 % Cumulative Flow Vol Bin
        if (MinStableTime < 0)
            [~,BestVolStep_ASCParamRow] = sortrows(GroupID_AllVidData_gCell{ugCt,5},[8 9],{'descend','ascend'}); 
        elseif (MinStableTime > -eps)
            [~,BestVolStep_ASCParamRow] = sortrows(GroupID_AllVidData_gCell{ugCt,5},[8 9],{'ascend','descend'}); 
        end
    end

    % Bypass Rows w/ Less than BinCtLim Complete Flow Vol Bins 
    BestVolStep_ASCParamRow(ismember(BestVolStep_ASCParamRow,find(GroupID_AllVidData_gCell{ugCt,5}(:,3)<=BinCtLlim)),:) = [];
    BestVolStep_ASCParamRow(ismember(BestVolStep_ASCParamRow,find(GroupID_AllVidData_gCell{ugCt,5}(:,8)<=MinStableVol)),:) = []; 
    
    if isempty(BestVolStep_ASCParamRow)
        continue;
    end

    BestVolStep_ASCParamRow = BestVolStep_ASCParamRow(1);
    BestVolStep = GroupID_AllVidData_gCell{ugCt,5}(BestVolStep_ASCParamRow,1);
    
    % Access RData corresponding to BestStepRow
    BestVolStep_RData = GroupID_AllVidData_gCell{ugCt,4}{1,1,BestVolStep_ASCParamRow};
    
    % Cum Flow Vol (nL) at which stable WBC conc begins for GroupID
    if UI.ASCOption==2  
        BestVolStep_CumVolStart = ...
            GroupID_AllVidData_gCell{ugCt,5}(BestVolStep_ASCParamRow,6) + GroupID_AllVidData_gCell{ugCt,5}(BestVolStep_ASCParamRow,1);
    end

    % Cum Imaging Time (s) at which stable WBC conc begins for GroupID
    BestVolStep_CumTimeStart = BestVolStep_RData(round(BestVolStep_RData(:,2),4)==round(BestVolStep_CumVolStart,4),4); 
    if (MinStableTime > -eps) & (BestVolStep_CumTimeStart < MinStableTime) 
        continue;
    end

    % Flow Vol Bin Idx at which stable WBC conc begins for GroupID
    BestVolStep_BinIdxStart = GroupID_AllVidData_gCell{ugCt,5}(BestVolStep_ASCParamRow,4);
    BestVolStep_BinIdxEnd = GroupID_AllVidData_gCell{ugCt,5}(BestVolStep_ASCParamRow,5);

    % Max. WBC Conc Ratio Difference from 1 among Flow Vol Bins during Stable WBC Conc
    BestVolStep_MaxRatioDiff = GroupID_AllVidData_gCell{ugCt,5}(BestVolStep_ASCParamRow,9);

    % Final Volume Bin's WBC Concentration
    BestVolStep_FinalVolBinConc = BestVolStep_RData(end,7);

    % Record BestVolStep_ASCParamRow in ASCParam
    GroupID_AllVidData_gCell{ugCt,5}(BestVolStep_ASCParamRow,10) = 1; 

    % Record WBC Stable Imaging Time and Cum Vol Range Data
    GroupID_AllVidData_gCell{ugCt,7} = BestVolStep_CumTimeStart;
    GroupID_AllVidData_gCell{ugCt,8} = BestVolStep_CumVolStart; 
    GroupID_AllVidData_gCell{ugCt,9} = BestVolStep; 
    GroupID_AllVidData_gCell{ugCt,10} = BestVolStep_BinIdxStart; 
    GroupID_AllVidData_gCell{ugCt,11} = BestVolStep_BinIdxEnd; 
    GroupID_AllVidData_gCell{ugCt,12} = BestVolStep_MaxRatioDiff; 
    GroupID_AllVidData_gCell{ugCt,13} = BestVolStep_FinalVolBinConc;  

end

clearvars BestVolStep*;
clearvars BinCtLlim ugCt;


%% Generate Vol Bin Concentration plots relevant to Auto Stability Criteria (ASC)

if UI.GroupIDPlotOption == 'y'

    [Plot_fPlotASCVolBinConc_gCell] = fPlotASCVolBinConc(GroupID_AllVidData_gCell, UI);

    for ugCt = 1:height(GroupID_AllVidData_gCell)
        if ~isempty(Plot_fPlotASCVolBinConc_gCell{ugCt,1})
            imwrite(Plot_fPlotASCVolBinConc_gCell{ugCt,1},strcat(SaveVBinPlotFilePath,'Plot_FlowVolBin_OB_StableWBCConcAnalysis_', GroupID_AllVidData_gCell{ugCt,1},'_',timestamp,'.tif'));
        end
    end

end

clearvars Plot_fPlotASCVolBinConc_gCell ugCt;


%% Convert Summary WBC Stability Start Time and Flow Vol Info to Export-able Table format; Save

fprintf('Converting Summary WBC Stability Start Time and Flow Vol Info to Export-able Table format; Saving ...\n');

GroupID_AllVidData_WBCConcStability_Tbl = GroupID_AllVidData_gCell(:,[1 3 7 8 9 12 13 14]); 
if UI.ASCOption == 2
    GroupID_AllVidData_WBCConcStability_Tbl = cell2table(GroupID_AllVidData_WBCConcStability_Tbl,...
        'VariableNames',{'GroupID','WinVesselDiameter_um','ASCImgTimeStart_s','ASCCumFlowVolStart_nL','FlowVolBinSize_nL','MaxWBCConcRatioDiff','FinalCumVolBinWBCConc_KuL','REFWBCConc_KuL'});
end

GroupID_AllVidData_WBCConcStability_Tbl.ASCOption = repmat(UI.ASCOption,height(GroupID_AllVidData_WBCConcStability_Tbl),1); 

GroupID_AllVidData_WBCConcStability_Tbl.TotalVolOB_nL = cellfun(@(x) max(x(:,31),[],'all'),GroupID_AllVidData_gCell(:,2),'UniformOutput',true);
GroupID_AllVidData_WBCConcStability_Tbl.TotalVolHi_nL = cellfun(@(x) max(x(:,31:36),[],'all'),GroupID_AllVidData_gCell(:,2),'UniformOutput',true);

GroupID_AllVidData_WBCConcStability_Tbl.TotalImgTime_s = cellfun(@(x) max(x(:,3),[],'all'),GroupID_AllVidData_gCell(:,2),'UniformOutput',true);

ZScr_50prc_Temp = cell2mat(cellfun(@(x) median(cell2mat(x),1),GroupID_AllVidData_gCell(:,6),'UniformOutput',false)); 
GroupID_AllVidData_WBCConcStability_Tbl.WBCZScr50prc = ZScr_50prc_Temp(:,2); 
ZScr_75prc_Temp = cell2mat(cellfun(@(x) prctile(cell2mat(x),75,1),GroupID_AllVidData_gCell(:,6),'UniformOutput',false)); 
GroupID_AllVidData_WBCConcStability_Tbl.WBCZScr75prc = ZScr_75prc_Temp(:,2); 

GroupID_AllVidData_WBCConcStability_Tbl.MeanVelocityOB_ums = cellfun(@(x) mean(x(:,43),'all'),GroupID_AllVidData_gCell(:,2),'UniformOutput',true);

writetable(GroupID_AllVidData_WBCConcStability_Tbl,strcat(SaveFilePath,'GroupID_AllVidData_WBCConcStability_Tbl_',timestamp,'.csv'));

clearvars x_plot y_plot ZScr_*prc_Temp; 


%% Plot Statistics of Starting ImgTime (s) and CumFlowVol (nL) at which Stable WBC Conc begins.

if any(GroupID_AllVidData_WBCConcStability_Tbl.ASCImgTimeStart_s>eps) | any(GroupID_AllVidData_WBCConcStability_Tbl.ASCCumFlowVolStart_nL>eps)
    [Plot_ASCVolBinStats] = fPlotASCVolBinStats(GroupID_AllVidData_WBCConcStability_Tbl,UI);
    imwrite(Plot_ASCVolBinStats,strcat(SaveFilePath,'Plot_BoxPlot_StableWBCConc_StartImgTimeCumVol_',timestamp,'.tif'));
end

clearvars Plot_ASCVolBinStats;


%% Remove GroupIDs not meeting ASC from GroupID_CumTimevsCumConc_PlotData_tlCell

[GroupID_CumTimevsCumConc_PlotData_tlCell_mCell,...
    Tbl_GroupID_CumTimevsCumConc_Stable_lmCell, TblLabel_GroupID_CumTimevsCumConc_Stable_lmCell] = ...
    fCleanASCPlotData(GroupID_CumTimevsCumConc_PlotData_tlCell_mCell,GroupID_AllVidData_WBCConcStability_Tbl,UI,1);

% Save table
% % % for torvlCt = 1:height(Tbl_GroupID_CumTimevsCumConc_Stable_lmCell) 
% % %     for mCt = 1:1
% % %         writetable(Tbl_GroupID_CumTimevsCumConc_Stable_lmCell{torvlCt,1,mCt},strcat(SaveTBinAllPlotCumTimevsCumConcFilePath,'Data_CumTimevsCumWBCConc_ASC',TblLabel_GroupID_CumTimevsCumConc_Stable_lmCell{torvlCt,1,mCt},'s_',timestamp,'.csv'))
% % %     end
% % % end

clearvars Tbl_GroupID_CumTimevsCumConc_Stable_lmCell TblLabel_GroupID_CumTimevsCumConc_Stable_lmCell;


%% Remove GroupIDs not meeting ASC from GroupID_CumVolvsCumConc_PlotData_vlCell

[GroupID_CumVolvsCumConc_PlotData_vlCell_mCell,...
    Tbl_GroupID_CumVolvsCumConc_Stable_lmCell, TblLabel_GroupID_CumVolvsCumConc_Stable_lmCell] = ...
    fCleanASCPlotData(GroupID_CumVolvsCumConc_PlotData_vlCell_mCell,GroupID_AllVidData_WBCConcStability_Tbl,UI,2);

% Save table
% % % for torvlCt = 1:height(Tbl_GroupID_CumVolvsCumConc_Stable_lmCell) % # Cum time limit
% % %     for mCt = 1:1 
% % %         writetable(Tbl_GroupID_CumVolvsCumConc_Stable_lmCell{torvlCt,1,mCt},strcat(SaveTBinAllPlotCumVolvsCumConcFilePath,'Data_CumVolvsCumWBCConc_ASC',TblLabel_GroupID_CumVolvsCumConc_Stable_lmCell{torvlCt,1,mCt},'nL_',timestamp,'.csv'))
% % %     end
% % % end

clearvars Tbl_GroupID_CumVolvsCumConc_Stable_lmCell TblLabel_GroupID_CumVolvsCumConc_Stable_lmCell;


%% REPEAT Plot: Cum Imaging Time (sec) (x) vs WBC Concentration (K/uL) (y) for each GroupID (Vessel Identifier)

% % % if any(GroupID_AllVidData_WBCConcStability_Tbl.ASCImgTimeStart_s>eps)
% % % 
% % %     [Plot_GroupID_CumTimevsCumConc_lpmCell, PlotLabel_GroupID_CumTimevsCumConc_lpmCell] = ...
% % %         fPlotCumWBCConcASC(GroupID_CumTimevsCumConc_PlotData_tlCell_mCell,UI,1); 
% % % 
% % %     % Save plot
% % %     for torvlCt = 1:height(Plot_GroupID_CumTimevsCumConc_lpmCell) 
% % %         for pCt = 1:1
% % %             for mCt = 1:1
% % %                 imwrite(Plot_GroupID_CumTimevsCumConc_lpmCell{torvlCt,pCt,mCt},strcat(SaveTBinAllPlotCumTimevsCumConcFilePath,'Plot_CumTimevsCumWBCConc_ASC',PlotLabel_GroupID_CumTimevsCumConc_lpmCell{torvlCt,pCt,mCt},'s_',timestamp,'.tif'));
% % %             end
% % %         end
% % %     end
% % % 
% % % end

clearvars Plot*_GroupID_CumTimevsCumConc_lpmCell torvlCt pCt mCt;


%% REPEAT Plot: Cum Blood Vol (nL) (x) vs Cum WBC Concentration (K/uL) (y) for each GroupID (Vessel Identifier)

% % % if any(GroupID_AllVidData_WBCConcStability_Tbl.ASCCumFlowVolStart_nL>eps)
% % % 
% % %     % = Plot
% % %     [Plot_GroupID_CumVolvsCumConc_lpmCell, PlotLabel_GroupID_CumVolvsCumConc_lpmCell] = ...
% % %         fPlotCumWBCConcASC(GroupID_CumVolvsCumConc_PlotData_vlCell_mCell,UI,2); 
% % % 
% % %     % Save plot
% % %     for torvlCt = 1:height(Plot_GroupID_CumVolvsCumConc_lpmCell) % # Cum volume limit
% % %         for pCt = 1:1
% % %             for mCt = 1:1
% % %                 imwrite(Plot_GroupID_CumVolvsCumConc_lpmCell{torvlCt,pCt,mCt},strcat(SaveTBinAllPlotCumVolvsCumConcFilePath,'Plot_CumVolvsCumWBCConc_ASC',PlotLabel_GroupID_CumVolvsCumConc_lpmCell{torvlCt,pCt,mCt},'nL_',timestamp,'.tif'));
% % %             end
% % %         end
% % %     end
% % % 
% % % end

clearvars Plot*_GroupID_CumVolvsCumConc_lpmCell torvlCt pCt mCt;


%% Generate WBC Concentration Ratio (Final Cumulative Vol Bin Conc / REF Conc) boxplot

[Plot_ConcRatio_Cell] = fPlotASCConcRatio(GroupID_AllVidData_WBCConcStability_Tbl,UI);
imwrite(Plot_ConcRatio_Cell{:,:,1},strcat(SaveFilePath,'Plot_BoxPlot_WBCConcRatio_',timestamp,'.tif'));


%% Generate Cell Ct Window Vessel Diameter Histogram

[Plot_WinSz] = fPlotWinVesselDiameterDistrb(GroupID_AllVidData_WBCConcStability_Tbl,UI);
imwrite(Plot_WinSz,strcat(SaveFilePath,'Plot_HistCDF_WinVesselDiameter_',timestamp,'.tif'));


%% Organize Output

oASC.GroupID_CumTimevsCumConc_mCell = GroupID_CumTimevsCumConc_mCell; 
oASC.GroupID_CumVolvsCumConc_mCell = GroupID_CumVolvsCumConc_mCell; 
oASC.GroupID_CumVolvsCumTime_mCell = GroupID_CumVolvsCumTime_mCell;  

% Summary Table
oASC.GroupID_AllVidData_WBCConcStability_Tbl = GroupID_AllVidData_WBCConcStability_Tbl;

oASC.GroupID_CumTimevsCumConc_PlotData_tlCell_mCell = GroupID_CumTimevsCumConc_PlotData_tlCell_mCell;
oASC.GroupID_CumVolvsCumConc_PlotData_vlCell_mCell = GroupID_CumVolvsCumConc_PlotData_vlCell_mCell;

oASC.GroupID_AllVidData_gCell = GroupID_AllVidData_gCell; 
oASC.GroupID_Data_fCell = GroupID_Data_fCell; 

clearvars GroupID_CumTimevsCumConc_mCell GroupID_CumVolvsCumConc_mCell GroupID_CumVolvsCumTime_mCell;
clearvars GroupID_AllVidData_WBCConcStability_Tbl;
clearvars GroupID_CumTimevsCumConc_PlotData_tlCell_mCell GroupID_CumVolvsCumConc_PlotData_vlCell_mCell;
clearvars GroupID_AllVidData_gCell GroupID_Data_fCell;
clearvars GroupID_VctyStat;


%% Delete parameters unlikely to be useful after this module

clearvars Plot_ConcRatio Plot_WinSz;
clearvars Save*FilePath -except SaveFilePath;


 %% Save Workspace

 fprintf('Saving Workspace...\n');
 save(strcat(SaveFilePath,'Workspace_WBCConcStabilityAnalysis_',timestamp,'.mat'),'-v7.3');
 

%% Save script in directory
% html format; do not evaulate code or save figures

ScriptName=mfilename;
PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
publish(strcat(ScriptName,'.m'),PublishOptions);

