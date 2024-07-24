clc; close all hidden; fclose('all'); clearvars;

BatchRun.Var00_FoldernameString = '/Users/j/Documents/MATLAB/LeukoCt/Data_In/VesselROI_In/'; % End w/ "/"
BatchRun.OutputFoldernameString = '/Users/j/Documents/MATLAB/LeukoCt/Data_Out/Pipeline_BATCHAUTO/'; % End w/ "/"

BatchRun.Script_List = ...
    {'STEP01_VesselROI.m',...
    'STEP02_SCS_VesselSz.m',...
    'STEP03_Velocity.m',...
    'STEP04_Velocity_Volume.m',...
    'STEP05_WBCCount.m',...
    'STEP06_WBCConc.m',...
    'STEP07_Report.m'};


%% Prepare variables relevant to all files

Var00_FolderFileList = dir(fullfile(BatchRun.Var00_FoldernameString,'*.mat'));

% = Prepare Summary Cell
% Col 1: Sample ID
% Col 2: 1 if Successful Run for Script in Script_List (Col# in Array), 0 otherwise
% Col 3: timestamp of current run for SampleID
% Col 4: timestamp of manual ROI draw for SampleID

BatchRunSummary_fCell = ... 
    repmat({"",zeros(1,numel(BatchRun.Script_List),'logical'),"",""},height(Var00_FolderFileList),1);


%%

for fCt = 1:height(Var00_FolderFileList)
    
    clearvars -except BatchRunData_Cell_All r_Data_Label_Cell BatchRunSummary_fCell BatchRun Var00_FolderFileList fCt;

    %% Run MOD00: Vessel ROI
    % New timestamp is given to Auto Run

    % = Extract SampleIDString
    expr_SampleID = '(?-i)\d+_P\d+_V\d+_\d+(_S\d+-\d+)?(_[a-z])?'; % Sample ID, Optional Slice & vessel, ex 20231005_P15_V01_01_S01-99_a
    [StartIdx,EndIdx] = regexp(Var00_FolderFileList(fCt).name,expr_SampleID);
    BatchRunSummary_fCell{fCt,1} = Var00_FolderFileList(fCt).name(StartIdx:EndIdx); % Col1: SampleID
    clearvars StartIdx EndIdx expr_SampleID;
    
    % = Extract OBM Video Input Filename
    Var00_FilenameString = fullfile(Var00_FolderFileList(fCt).folder,Var00_FolderFileList(fCt).name);
    load(Var00_FilenameString,'ImgStack_H5orTIF_FilenameString');

    try
        MOD00_AutoRun_Option = 1;
        run(BatchRun.Script_List{1});
        BatchRunSummary_fCell{fCt,3} = timestamp;
    catch
        continue;
    end
    
    % = Extract manual ROIDraw (MOD00) related variables to load within saved
    % MAT file; Load
    Var00 = who('-file',Var00_FilenameString); 
    Var00_LoadStr = ("PxLength"|"Crop"|"PolygonROI_startXY"|"TimeParam"|"tFr_MeanTissueIntensity"|"oROI");
    Var00_Load = Var00(cellfun(@(x) ~isempty(strfind(x,Var00_LoadStr)),Var00,'UniformOutput',true));
    load(Var00_FilenameString,Var00_Load{:}); 
    ROIManDraw = load(Var00_FilenameString,'timestamp'); 
    BatchRunSummary_fCell{fCt,4} = ROIManDraw.timestamp;
    clearvars Var00_FilenameString Var00 Var00_LoadStr Var00_Load ROIManDraw;

    % = Save STEP01 (AUTO) Workspace
    fprintf('Saving STEP01 (AUTO) Workspace...\n');
    save(strcat(SaveMODFilePath,'Workspace_STEP01Var_',SampleIDString,'_AUTO_',timestamp,'.mat'),'-regexp','^(?!(fCt)$).');
    varList_MOD00 = who;
    varList_All = varList_MOD00;

    % = Mark successful run 
    BatchRunSummary_fCell{fCt,2}(1,1) = 1; 

    clearvars MOD00_AutoRun_Option; % Important


    %% Run rest of pipeline

    for sCt = 2:numel(BatchRun.Script_List) 
        try
            run(BatchRun.Script_List{sCt});
            BatchRunSummary_fCell{fCt,2}(1,sCt) = 1; 
        catch
            break; 
        end
    end

    clearvars sCt;


    %% Append Report_Cell_All

    if ~exist("BatchRunData_Cell_All") & (sum(BatchRunSummary_fCell{fCt,2})==numel(BatchRun.Script_List))
        BatchRunData_Cell_All = r_Data_Cell; 
    elseif exist("BatchRunData_Cell_All") & (sum(BatchRunSummary_fCell{fCt,2})==numel(BatchRun.Script_List))
        BatchRunData_Cell_All = vertcat(BatchRunData_Cell_All,r_Data_Cell);
    elseif (sum(BatchRunSummary_fCell{fCt,2})<numel(BatchRun.Script_List)) 
        rmdir(SaveFilePath,'s'); 
    end

    clearvars r_Data_Cell;

end


%% Save

% = Prepare Save Folder; Summary Files save one subfolder layer above
% (outside sample folder)
slashIdx   = strfind(SaveFilePath,'/');
SaveFilePathPrev = convertStringsToChars(SaveFilePath);
SaveFilePathPrev = strcat(convertCharsToStrings(SaveFilePathPrev(1:slashIdx(end-1)-1)),'/');

writecell(BatchRunSummary_fCell,strcat(SaveFilePathPrev,'BatchRun_Summary_',BatchRunSummary_fCell{1,3},'-',BatchRunSummary_fCell{end,3},'.csv'));

if exist("BatchRunData_Cell_All")
    writecell(BatchRunData_Cell_All,strcat(SaveFilePathPrev,'BatchRun_Data_',BatchRunSummary_fCell{1,3},'-',BatchRunSummary_fCell{end,3},'.csv'));
    writecell(vertcat(r_Data_Label_Cell,BatchRunData_Cell_All),strcat(SaveFilePathPrev,'BatchRun_Data_wLabel_',BatchRunSummary_fCell{1,3},'-',BatchRunSummary_fCell{end,3},'.csv'));
end

clearvars slashIdx;

