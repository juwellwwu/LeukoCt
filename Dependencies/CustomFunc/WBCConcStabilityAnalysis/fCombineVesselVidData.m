
% fCombineVesselVidData() combines cell ct pipeline data of all videos
% of one vessel. "vidDirLlist_In" is a dir() struct object (w/ fields =
% name, folder, data, bytes, isdir, datenum). 
% GateXY_Method_Report_In is the cell count method used (= 4 "PxWin")

% "tMidPtData_Out" has following columns:

% 1) VidIdx: Which video

% 2) tMidPt_FrStart (s, referenced to start of video 1)
% 3) tMidPt_FrEnd (s, referenced to start of video 1) 
% 4) tMidPt_FrLength (s, referenced to start of video 1) 
% 5) tMidPt_FrCtr (s, referenced to start of video 1)

% 6) Incremental Blood Flow Vol (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 7) Incremental Blood Flow Vol StdEr (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 8) Incremental Blood Flow Vol (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 9) Incremental Blood Flow Vol StdEr (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK 

% 10) WBC Count between tMidPt_FrStart and tMidPt_FrEnd 
% 11) WBC Z-Score between tMidPt_FrStart and tMidPt_FrEnd 

% 12) Incremental CellConc (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 13) Incremental CellConc Lo (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 14) Incremental CellConc Hi (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 

% 15) Incremental CellConc (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 16) Incremental CellConc Lo (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 17) Incremental CellConc Hi (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 

% 18) XYWin CellCt Window Diameter between tMidPt_FrStart & tMidPt_FrEnd 

% 19) Flag 01: Radius [DUMMY]
% 20) Flag 02: VctyEr, OB
% 21) Flag (Metric) 02: VctyEr, OB 
% 22) Flag 03: VctyEr, SK 
% 23) Flag (Metric) 03: VctyEr, SK 
% 24) Flag 04: FlowStb [DUMMY]
% 25) Flag 05: RBCBkgd [DUMMY]
% 26) Flag 06: CellCt, by GroupID [DUMMY]
% 27) Flag (Metric) 06: CellCt, by GroupID [DUMMY]

% 28) Flag X: for future use
% 29) Flag X: by future use

% 30) Flag ALL: which tMidPt to cut considering all flags above

% 31) Cumulative Blood Flow Vol (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 32) Cumulative Blood Flow Vol Hi (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 33) Cumulative Blood Flow Vol Lo (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB

% 34) Cumulative Blood Flow Vol (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 35) Cumulative Blood Flow Vol Hi (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 36) Cumulative Blood Flow Vol Lo (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK 

% 37) Cumulative CellConc (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 38) Cumulative CellConc Lo (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 39) Cumulative CellConc Hi (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 

% 40) Cumulative CellConc (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 41) Cumulative CellConc Lo (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 42) Cumulative CellConc Hi (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 

% 43) Incremental Blood Flow Velocity (um/s) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 44) Incremental Blood Flow Velocity StdEr (um/s) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 45) Incremental Blood Flow Velocity (um/s) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 46) Incremental Blood Flow Velocity StdEr (um/s) between tMidPt_FrStart & tMidPt_FrEnd, SK 

% 47) Vessel Lumen Diameter (um) between tMidPt_FrStart & tMidPt_FrEnd 
% 48) Vessel Flow Profile Diameter (um) between tMidPt_FrStart & tMidPt_FrEnd 

% 49) Expected False+ WBC Count from RBC Bkgd between tMidPt_FrStart and tMidPt_FrEnd 
% 50) Corrected WBC Count between tMidPt_FrStart and tMidPt_FrEnd 

% "CellLocZScore_Out" has the information of WBC locations, Z-Score
% Col 1: WBC Loc, referenced to start of all videos after cut
% Col 2: WBC Z Score

% "RefConc_KuL_Out"  = CBC WBC Concentration, in K/uL

function [tMidPtData_Out,CellLocZScore_Cell_Out,RefConc_KuL_Out] = fCombineVesselVidData(vidDirStruct_In,GateXY_Method_Report_In,UI_In)


FlagCutOption_In = UI_In.FlagtMidPtCutOption;
FlagCellCt_Diameter_um_In = UI_In.FlagtMidPtCut.CellCtDiameterLim_um;

CtCorrOption_In = UI_In.CtCorrOption;


%%

RowCt = 0; 
tMidPtData_Out = -9999*ones(5000,49,'single'); 

CellLocZScore_Cell_Out = cell(0,2); 

for vCt = 1:height(vidDirStruct_In)

    % = Load necessary variables from MAT file
    matFileList = dir(fullfile(vidDirStruct_In(vCt).folder,vidDirStruct_In(vCt).name,'*/*.mat'));
    for mfCt = 1:height(matFileList) 
        mat_FilenameString = fullfile(matFileList(mfCt).folder,matFileList(mfCt).name);
        mat_Var = who('-file',mat_FilenameString);
        mat_Var_LoadStr = ... 
            ("PxLength" | "TimeParam" | "oSCS" | "oRadius" | "oVelocity" | "oVolume" | "oCellCt" | "oCellConc" | "oFlag");
        mat_Var_Load = mat_Var(cellfun(@(x) ~isempty(strfind(x,mat_Var_LoadStr)),mat_Var,'UniformOutput',true));
        if ~isempty(mat_Var_Load) 
            load(mat_FilenameString,mat_Var_Load{:}); 
        end
    end
    
    RefConc_KuL_Out = oCellConc_mCell{1,4}.Ref_KuL;

    clearvars mat_FilenameString mat_Var mat_Var_LoadStr mat_Var_Load matFileList mfCt;

    % === 1) VidComboData_Out

    % = Which cell ct window and its vessel block to extract data from loaded outputs
    BrightBlkorWinIdx = find(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.BrightRank==1); 
    BlkSCS_BrightBlkorWinIdx = oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.WinSkelPxBlk(find(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.BrightRank==1)); 

    rS = RowCt+1; 
    rE = RowCt+numel(TimeParam.tMidPt_FrLength); 

    % = Load columns of Data_Out that requires combining video data

    % 1) VidIdx: Which video
    tMidPtData_Out(rS:rE,1) = repmat(vCt,numel(TimeParam.tMidPt_FrLength),1);

    % 2) tMidPt_FrStart (s, referenced to start of video 1)
    % 3) tMidPt_FrEnd (s, referenced to start of video 1) 
    % 4) tMidPt_FrLength (s, referenced to start of video 1) 
    % 5) tMidPt_FrCtr (s, referenced to start of video 1) 

    tMidPtData_Out(rS:rE,4) = TimeParam.tMidPt_FrLength.*Z_PxLength; % will update Col 2,3,5 later to change start reference

    % 6) Incremental Blood Flow Vol (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
    % 7) Incremental Blood Flow Vol StdEr (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB
    % 8) Incremental Blood Flow Vol (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK
    % 9) Incremental Blood Flow Vol StdEr (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK

    tMidPtData_Out(rS:rE,6) = cell2mat(oVolume.V.tMidPt_AllSkOBwMean_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength^3./1E6;
    tMidPtData_Out(rS:rE,7) = cell2mat(oVolume.V.tMidPt_AllSkOBwStdEr_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength^3./1E6;
    tMidPtData_Out(rS:rE,8) = cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkMean_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength^3./1E6;
    tMidPtData_Out(rS:rE,9) = cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkStdEr_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength^3./1E6;

   % 10) WBC Count between tMidPt_FrStart and tMidPt_FrEnd 
   tMidPtData_Out(rS:rE,10) = cell2mat(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx}); % this parameter has no SCS vessel blk info

    % 11) (Mean) WBC Z-Score between tMidPt_FrStart and tMidPt_FrEnd 
    tMidPtData_Out(rS:rE,11) = cell2mat(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.tMidPt_MeanZScore_Cell{1,BrightBlkorWinIdx}); % this parameter has no SCS vessel blk info
    
    % 12) Incremental CellConc (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB
    % 13) Incremental CellConc Lo (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
    % 14) Incremental CellConc Hi (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB

    tMidPtData_Out(rS:rE,12) = cell2mat(oCellConc_mCell{1,GateXY_Method_Report_In}.ic_KuL.tMidPt_AllSkOBwMean_Cell)';
    tMidPtData_Out(rS:rE,13) = cell2mat(oCellConc_mCell{1,GateXY_Method_Report_In}.ic_KuL.tMidPt_AllSkOBwMean_Lo_Cell)';
    tMidPtData_Out(rS:rE,14) = cell2mat(oCellConc_mCell{1,GateXY_Method_Report_In}.ic_KuL.tMidPt_AllSkOBwMean_Hi_Cell)';

    % 15) Incremental CellConc (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
    % 16) Incremental CellConc Lo (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
    % 17) Incremental CellConc Hi (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 

    tMidPtData_Out(rS:rE,15) = cell2mat(oCellConc_mCell{1,GateXY_Method_Report_In}.ic_KuL.tMidPt_1SkOBwMean_SkMean_Cell)';
    tMidPtData_Out(rS:rE,16) = cell2mat(oCellConc_mCell{1,GateXY_Method_Report_In}.ic_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell)';
    tMidPtData_Out(rS:rE,17) = cell2mat(oCellConc_mCell{1,GateXY_Method_Report_In}.ic_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell)';

    % 18) XYWin CellCt Window Diameter (um) between tMidPt_FrStart & tMidPt_FrEnd 
    tMidPtData_Out(rS:rE,18) = repmat(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.WinVesselDiameter_px(1,BrightBlkorWinIdx).*XY_PxLength,...
        numel(TimeParam.tMidPt_FrLength),1); % may have different value for different videos due to different window location 

    % 19) Flag 01: Radius (DUMMY)
    % 20) Flag 02: VctyEr, OB 
    % 21) Flag (Metric) 02: VctyEr, OB 
    % 22) Flag 03: VctyEr, SK 
    % 23) Flag (Metric) 03: VctyEr, SK 
    % 24) Flag 04: FlowStb (DUMMY)
    % 25) Flag 05: RBCBkgd (DUMMY)
    % 26) Flag 06: CellCt, by GroupID (DUMMY)
    % 27) Flag (Metric) 06, by GroupID (DUMMY)

    tMidPtData_Out(rS:rE,19) = cell2mat(oFlag_Radius.L.tMidPt_1BlkAllSk_Flag_Cell(1,:,BlkSCS_BrightBlkorWinIdx))';

    tMidPtData_Out(rS:rE,20) = cell2mat(oFlag_VctyEr.OB.tMidPt_1BlkAllSk_Flag_Cell(1,:,BlkSCS_BrightBlkorWinIdx))';
    tMidPtData_Out(rS:rE,21) = cell2mat(oFlag_VctyEr.OB.tMidPt_1BlkAllSk_Metric_Cell(1,:,BlkSCS_BrightBlkorWinIdx))';

    tMidPtData_Out(rS:rE,22) = cell2mat(oFlag_VctyEr.SK.tMidPt_1BlkAllSk_Flag_Cell(1,:,BlkSCS_BrightBlkorWinIdx))';
    tMidPtData_Out(rS:rE,23) = cell2mat(oFlag_VctyEr.SK.tMidPt_1BlkAllSk_Metric_Cell(1,:,BlkSCS_BrightBlkorWinIdx))';

    tMidPtData_Out(rS:rE,24) = cell2mat(oFlag_FlowStb.tMidPt_1BlkAllSk_Flag_Cell(1,:,BlkSCS_BrightBlkorWinIdx))';

    tMidPtData_Out(rS:rE,25) = cell2mat(oFlag_RBCBkgd_mCell{1,GateXY_Method_Report_In}.tMidPt_1BlkAllSk_Flag_Cell(1,:,1))'; 

    tMidPtData_Out(rS:rE,26) = cell2mat(oFlag_CellCt_mCell{1,GateXY_Method_Report_In}.tMidPt_1BlkAllSk_Flag_Cell(1,:,BlkSCS_BrightBlkorWinIdx))'; 
    tMidPtData_Out(rS:rE,27) = cell2mat(oFlag_CellCt_mCell{1,GateXY_Method_Report_In}.tMidPt_1BlkAllSk_Metric_Cell(1,:,BlkSCS_BrightBlkorWinIdx))'; 
   
    % 28) Flag X: for future use
    % 29) Flag X: by future use

    % Velocity (um/s) Data
    % 43) Incremental Blood Flow Velocity (um/s) between tMidPt_FrStart & tMidPt_FrEnd, OB 
    % 44) Incremental Blood Flow Velocity StdEr (um/s) between tMidPt_FrStart & tMidPt_FrEnd, OB 
    % 45) Incremental Blood Flow Velocity (um/s) between tMidPt_FrStart & tMidPt_FrEnd, SK 
    % 46) Incremental Blood Flow Velocity StdEr (um/s) between tMidPt_FrStart & tMidPt_FrEnd, SK 

    tMidPtData_Out(rS:rE,43) = cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength./Z_PxLength;
    tMidPtData_Out(rS:rE,44) = cell2mat(oVelocity.tMidPt_AllSkOBwStdEr_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength./Z_PxLength;
    tMidPtData_Out(rS:rE,45) = cell2mat(oVelocity.tMidPt_1SkOBwMean_SkMean_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength./Z_PxLength;
    tMidPtData_Out(rS:rE,46) = cell2mat(oVelocity.tMidPt_1SkOBwMean_SkStdEr_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength./Z_PxLength;

    % Vessel Lumen and Flow Profile Diameter (um) Data
    % 47) Vessel Lumen Diameter (um) between tMidPt_FrStart & tMidPt_FrEnd 
    % 48) Vessel Flow Profile Diameter (um) between tMidPt_FrStart & tMidPt_FrEnd
    tMidPtData_Out(rS:rE,47) = cell2mat(oRadius.L.tMidPt_AllSk_Diameter_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength; 
    tMidPtData_Out(rS:rE,48) = cell2mat(oRadius.F.tMidPt_AllSk_Diameter_Cell(1,:,BlkSCS_BrightBlkorWinIdx)').*XY_PxLength; 

    % False-positive corrected WBC counts
    % 49) Expected False+ WBC Ct from RBC Bkgd
    % 50) Corrected WBC Count

   tMidPtData_Out(rS:rE,49) = ...
       oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.CellCtExpEr(1,BrightBlkorWinIdx) .* sum(tMidPtData_Out(rS:rE,10),"all") .* ...
       (TimeParam.tMidPt_FrLength./sum(TimeParam.tMidPt_FrLength,'all'))';

   if numel(find(isnan(tMidPtData_Out(rS:rE,49))))>0  
       if isfield(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt,'multScMAD_CellCt') 
           tMidPtData_Out(rS:rE,49) = (TimeParam.tMidPt_FrLength.*(1-normcdf(multScMAD_CellCt_In,0,1)))';
       else
            tMidPtData_Out(rS:rE,49) = (TimeParam.tMidPt_FrLength.*(1-normcdf(3.5,0,1)))';
       end
   end

   if CtCorrOption_In == 1
       tMidPtData_Out(rS:rE,50) =  tMidPtData_Out(rS:rE,10)-tMidPtData_Out(rS:rE,49); 
   elseif CtCorrOption_In == 0
        tMidPtData_Out(rS:rE,50) =  tMidPtData_Out(rS:rE,10);
   end

    % === 2) CellLocZScore_Cell_Out
    RelLoc = (oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.pkLoc{1,BrightBlkorWinIdx} - ...
        (TimeParam.tMidPt_FrStart(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.pktMidPt{1,BrightBlkorWinIdx}))' + 1) .* Z_PxLength;
    CellLocZScore_Cell_Out = vertcat(CellLocZScore_Cell_Out,...
        horzcat(mat2cell(RelLoc,cell2mat(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx})),...
        mat2cell(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.pkZScore{1,BrightBlkorWinIdx},cell2mat(oCellCt.Ct_mCell{1,GateXY_Method_Report_In}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx}))));

    % === Update vRowCt and ctMidPtCt
    RowCt = RowCt + numel(TimeParam.tMidPt_FrLength);
   
end

% Remove rows at which VidIdx (Col1) < 0, ie, not loaded with any data
tMidPtData_Out(tMidPtData_Out(:,1)<0,:) = [];

clearvars rS rE RowCt;


%% Option to Cut; Final Flag (combining all flags) in Col 30

% === Cut VidComboData_Out

% 30) Flag ALL: which tMidPt to cut considering all flags above
% VctyEr Flags
if FlagCutOption_In == 0 % Do not remove Flagged tMidPt of any flags
    tMidPtData_Out(:,30) = repmat(single(0),height(tMidPtData_Out),1);
elseif FlagCutOption_In == 1 % Remove tMidPt flagged by VctyEr_OB only
    tMidPtData_Out(:,30) = tMidPtData_Out(:,20);
elseif FlagCutOption_In == 2 % Remove tMidPt flagged by VctyEr_SK only
    tMidPtData_Out(:,30) = tMidPtData_Out(:,22);
elseif (FlagCutOption_In == 3) | (FlagCutOption_In == 4) | (FlagCutOption_In == 5)  % Remove tMidPt flagged by VctyEr_OB & VctyEr_SK
    tMidPtData_Out(:,30) = single(logical(tMidPtData_Out(:,20)) | logical(tMidPtData_Out(:,22)));
end

% RBC Bkgd Flag [DUMMY]
if FlagCutOption_In == 5 
    tMidPtData_Out(:,30) = single(logical(tMidPtData_Out(:,30)) | logical(tMidPtData_Out(:,25)));
end

FlagID01 = logical(tMidPtData_Out(:,30)); % assign FlagID for CellLocZScore_Cell_Out cut later
tMidPtData_Out(FlagID01,:) = [];

% = Update CellCt Flag [DUMMY] from by video to by GroupID
% 26) Flag 06: CellCt, by GroupID  (Original NONE)
% 27) Flag (Metric) 06: CellCt, by GroupID (Original NONE)
tMidPtData_Out(:,27) = repmat(sum(tMidPtData_Out(:,49),'all')./sum(tMidPtData_Out(:,10),'all'),height(tMidPtData_Out),1); % same value for all tMidPt of GroupID
tMidPtData_Out(:,26) = 0; % Reset
if tMidPtData_Out(1,27) > 0.1
     tMidPtData_Out(:,26) = single(1);
end

% 30) Flag ALL: which tMidPt to cut considering all flags above
% Low CellCt Flag
if (FlagCutOption_In == 4) | (FlagCutOption_In == 5) 
    tMidPtData_Out(:,30) = single(logical(tMidPtData_Out(:,30)) | (logical(tMidPtData_Out(:,26)) & (tMidPtData_Out(:,18)<(FlagCellCt_Diameter_um_In+eps))));
end

FlagID02 = logical(tMidPtData_Out(:,30)); % assign FlagID for CellLocZScore_Cell_Out cut later
tMidPtData_Out(FlagID02,:) = [];

% = Update Time Data Columns (2,3,5) such that reference to start of video 1
% 2) tMidPt_FrStart (s, referenced to start of video 1)
% 3) tMidPt_FrEnd (s, referenced to start of video 1) (Original Col 3)
% 4) tMidPt_FrLength (s, referenced to start of video 1) (Original NONE)
% 5) tMidPt_FrCtr (s, referenced to start of video 1) (Original Col 8)

tMidPtData_Out(:,2) = cumsum(vertcat(0,tMidPtData_Out(1:(end-1),4)));
tMidPtData_Out(:,3) = cumsum(tMidPtData_Out(:,4));
tMidPtData_Out(:,5) = 0.5*(tMidPtData_Out(:,2)+tMidPtData_Out(:,3));


% === Cut CellLocZScore_Cell_Out
CellLocZScore_Cell_Out(FlagID01,:) = [];
CellLocZScore_Cell_Out(FlagID02,:) = [];

% Switch relative WBC loc values to absolute WBC loc values (s)
for ceCt = 1:height(CellLocZScore_Cell_Out)
    if ~isempty(CellLocZScore_Cell_Out{ceCt,1})
        CellLocZScore_Cell_Out{ceCt,1} = CellLocZScore_Cell_Out{ceCt,1}+tMidPtData_Out(ceCt,2);
    end
end


clearvars ceCt FlagID*


%% Calculate columns of VidComboData_Out that requires combinED and cut video data

% 31) Cumulative Blood Flow Vol (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB
% 32) Cumulative Blood Flow Vol Hi (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB
% 33) Cumulative Blood Flow Vol Lo (nL) between tMidPt_FrStart & tMidPt_FrEnd, OB

tMidPtData_Out(:,31) = cumsum(tMidPtData_Out(:,6));
tMidPtData_Out(:,32) = cumsum(tMidPtData_Out(:,6)) + cumsum(tMidPtData_Out(:,7));
tMidPtData_Out(:,33) = cumsum(tMidPtData_Out(:,6)) - cumsum(tMidPtData_Out(:,7));

% 34) Cumulative Blood Flow Vol (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK
% 35) Cumulative Blood Flow Vol Hi (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 36) Cumulative Blood Flow Vol Lo (nL) between tMidPt_FrStart & tMidPt_FrEnd, SK 

tMidPtData_Out(:,34) = cumsum(tMidPtData_Out(:,8));
tMidPtData_Out(:,35) = cumsum(tMidPtData_Out(:,8)) + cumsum(tMidPtData_Out(:,9));
tMidPtData_Out(:,36) = cumsum(tMidPtData_Out(:,8)) - cumsum(tMidPtData_Out(:,9));

% 37) Cumulative CellConc (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 38) Cumulative CellConc Lo (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 
% 39) Cumulative CellConc Hi (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, OB 

tMidPtData_Out(:,37) = cumsum(tMidPtData_Out(:,50)) ./ tMidPtData_Out(:,31);
tMidPtData_Out(:,38) = cumsum(tMidPtData_Out(:,50)) ./ tMidPtData_Out(:,32);
tMidPtData_Out(:,39) = cumsum(tMidPtData_Out(:,50)) ./ tMidPtData_Out(:,33);

% 40) Cumulative CellConc (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK 
% 41) Cumulative CellConc Lo (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK
% 42) Cumulative CellConc Hi (K/uL) between tMidPt_FrStart & tMidPt_FrEnd, SK

tMidPtData_Out(:,40) = cumsum(tMidPtData_Out(:,50)) ./ tMidPtData_Out(:,34);
tMidPtData_Out(:,41) = cumsum(tMidPtData_Out(:,50)) ./ tMidPtData_Out(:,35);
tMidPtData_Out(:,42) = cumsum(tMidPtData_Out(:,50)) ./ tMidPtData_Out(:,36);

