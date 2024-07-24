clc; close all hidden; fclose('all');


%% User-specified variables

% === Input Image Stack to be analyzed
if ~exist("MOD00_AutoRun_Option") | isempty(MOD00_AutoRun_Option) | ~MOD00_AutoRun_Option

    clearvars;
    MOD00_AutoRun_Option = 0;    
    ImgStack_H5orTIF_FilenameString = ... % User input here *****
        '/Users/j/Documents/MATLAB/LeukoCt/Data_In/REGVids_In/20230418_P05_V10_05_1.72pxum_FastRegXYZIntEq.tif';

    % === Option to append SampleIDString
    % Use if use same video for >1 ROI analysis (ex. analyse 2 vessels in 1 video)
    % AppendSampleIDStr.Str is appended to SampleIDStr
    % NOTE: If file has manual cell count info, copy manual cell count file and rename so that
    % file names matched appended file names. This allows manual count
    % for multiple different vessels.
    AppendSampleIDStr.Option = 'n';
    AppendSampleIDStr.Str = 'a';

    % === Output Root Folder (include / at end) *****
    Output_RootFoldernameString = '/Users/j/Documents/MATLAB/LeukoCt/Data_Out/Pipeline_MANUAL/';

else

    % === Output Root Folder (include / at end) *****
    Output_RootFoldernameString = '/Users/j/Documents/MATLAB/LeukoCt/Data_Out/Pipeline_BATCHAUTO/';

end

% === OBM Video Properties *****
% Excel worksheet containing  pixel dimensions XY_PxLength & Z_PxLength
VidProp_FilenameString = ...
    '/Users/j/Documents/MATLAB/LeukoCt/Data_In/OBMVideoProperties_In.xlsx';

% === Reference CBC WBC Concentration (K/uL) *****
CBC_FilenameString = ...
    '/Users/j/Documents/MATLAB/LeukoCt/Data_In/CBC_In.xlsx';

% === Manual Cell Count (OPTIONAL) *****
% List of frames with manually recorded WBC position
% Auto search for file with same SampleID in folder
Fr_ManualCellCt_FoldernameString = '/Users/j/Documents/MATLAB/LeukoCt/Data_In/ManualCellCt_In/';

% === Option to modify pixel length in XY (um)
AdjMagParam.Option = 'n'; % default: 'n'
AdjMagParam.Target_XY_PxLength=1/1.72;  

% === Option to crop ImgStack in Z dimension
ZCropParam.Option = 'n'; % if "y", crop ImgStack in Z (time) dimension
ZCropParam.ZStart = 1; % Starting slice of ZCrop; Slice included in post-crop stack
ZCropParam.ZEnd = 1600; % Ending slice of ZCrop; Slice included in post-crop stack

% === Frames segments parameters
tSeg_FrLength = 200; % default: 200

% === Vessel Block Length, in number of 5-um cell diameters
% Vessel block Length = NumCellDiameter_in_SkelBlockLength+/-ß2 
NumCellDiameter_in_SkelBlockLength = 10; % default: 10

% === Whether Vessel Diameter is treated as time dependent
Vessel_Diameter_TimeDependence = 'n'; % default: 'n'

% === Area of Cell Count Window
% Area expressed as fraction occupied by average WBC (12.0um) Area;
WBCtoGateXYWinAreaRatio_Init = [0.40 0.35 0.30 0.25 0.20];

% === Option to create DEMO / Illustration images
DEMOOption = 'n'; % default: 'n'


%% Load Video Parameters: XY_PxLength, Z_PxLength

fprintf('Loading pixel size information of ImgStack...\n');

expr_SampleID = '(?-i)\d+_P\d+_V\d+_\d+(?:_S\d+-\d+)?'; % Sample ID, Optional Slice, ex 20231005_P15_V01_01_S01-99
[SampleIDString, XY_PxLength, Z_PxLength, ~, ~, ~] = ...
    fExtractVidProperty(VidProp_FilenameString,ImgStack_H5orTIF_FilenameString, expr_SampleID);

if ~MOD00_AutoRun_Option
    if AppendSampleIDStr.Option == 'y'
        SampleIDString = horzcat(SampleIDString,'_',AppendSampleIDStr.Str);
    end
else
    SampleIDString = BatchRunSummary_fCell{fCt,1};
end

clearvars expr_SampleID;


%% Create Output Folder Names

% = Output folder (do NOT include / at end)
if MOD00_AutoRun_Option
    OutputFolder = strcat(Output_RootFoldernameString,SampleIDString,'_AUTO_PIPELINE');
else
    OutputFolder = strcat(Output_RootFoldernameString,SampleIDString,'_PIPELINE');
end

timestamp = string(datetime('now'),'yyMMddHHmm');
SaveFilePath = strcat(OutputFolder,'_',timestamp,'/');
mkdir(SaveFilePath);

SaveMODFilePath = strcat(SaveFilePath,'STEP01_VesselROI/');
mkdir(SaveMODFilePath);

if DEMOOption == 'y'
    SaveMODDEMOFilePath = strcat(SaveMODFilePath,'DEMOIllus/');
    mkdir(SaveMODDEMOFilePath);
end


%% ROI Draw

if ~MOD00_AutoRun_Option

    %% Load and prepare ImgStack for processing

    fprintf('Loading and preparing video for processing...\n');

    expr_SampleID = '(?-i)\d+_P\d+_V\d+_\d+(?:_S\d+-\d+)?'; % Sample ID, Optional Slice, ex 20231005_P15_V01_01_S01-99
    [ImgStack,XY_PxLength,ZCropParamMOD,SampleIDString] = ...
        fPrepTifVolume(ImgStack_H5orTIF_FilenameString,XY_PxLength,AdjMagParam,ZCropParam,tSeg_FrLength,expr_SampleID,'n');

    if AppendSampleIDStr.Option == 'y'
        SampleIDString = horzcat(SampleIDString, '_', AppendSampleIDStr.Str);
    end

    clearvars expr_SampleID SampleIDString2;


    %% Calculate Time Parameters

    TimeParam = fDefineTimeParamStruct(tSeg_FrLength,ZCropParamMOD);


    %%  Load Manual Cell Ct Frame information, if info exists in directory

    [Fr_ManualCellCt] = fLoadManualCellCtFr(Fr_ManualCellCt_FoldernameString,SampleIDString,ZCropParamMOD);


    %% Create preROIDraw Vessel Mask
    % * PreROI_PreClean: Vessel Lumen + Flow Profile Masks drawn by fCreateRegImgZThresh;
    % Vessel ROI not specified, which vessel to be used not specified, mask can
    % include unwanted areas.
    % * PreROI: Vessel ROI not specified. User specifies vessel to analyse with foreground (blue
    % line) draw. User can then use line eraser to remove unwanted areas.

    [Img_Zstdev,Img_Thresh_PreROI_PreClean_sCell,Img_Thresh_PreROI_sCell,Img_Skel_PreROI_sCell] = ...
        fCreatePreROIDrawVesselMask(ImgStack,XY_PxLength,SaveMODFilePath);


    %% Obtain approximate tissue intensity value

    TissueMask = (~max(cell2mat(permute(Img_Thresh_PreROI_PreClean_sCell,[3 2 1])),[],3)); 
    [tFr_MeanTissueIntensity] = fCalctFrTissueIntensity(ImgStack,TissueMask,'n');

    clearvars TissueMask*;


    %% Define ROI or Vessel Length to be analyzed

    fprintf('User-assisted ROI Draw...\n');

    % = Overlay ImgStack and Skeleton
    Img_ROIDraw = imfuse(0.50*(0.50*Img_Thresh_PreROI_sCell{1,1}+0.50*Img_Thresh_PreROI_sCell{2,1}),min(ImgStack,[],3),'Scaling','none','ColorChannels',[1 2 0]);
    Img_ROIDraw(find(Img_Skel_PreROI_sCell{1,1})) = 0.5;
    Img_ROIDraw(find(Img_Skel_PreROI_sCell{2,1})) = 0.5;

    FigHandle = figure;
    set(FigHandle, 'MenuBar', 'none');
    set(FigHandle, 'ToolBar', 'none');
    imshow(Img_ROIDraw,'InitialMagnification',300);

    % = User input: line position recorded as [x1 y1; x2 y2]
    fprintf('**USER INPUT**: draw ROI...\n');
    fprintf('Draw where skeleton looks correct. Initiate drawing near source of flow...\n');
    PolygonROI = drawpolygon('LineWidth',1,'Color','cyan');

    PolygonROI_startXY = PolygonROI.Position(1,:); % Keep 1st x,y position to record direction of flow

    % = Create ROI Mask; record where ROI is drawn
    Img_ROI = createMask(PolygonROI);
    imwrite(imfuse(bwperim(Img_ROI),Img_ROIDraw,'blend','Scaling','none'),strcat(SaveMODFilePath,'Img_PolygonROI_',timestamp,'.tif'));

    % = For Demo Illustrations only
    if DEMOOption == 'y'
        imwrite(uint8(255.*Img_Thresh_PreROI_sCell{1,1}),strcat(SaveMODDEMOFilePath,'Img_illus_lThresh_PreROI_',timestamp,'.tif'));
        imwrite(uint8(255.*bwperim(Img_Thresh_PreROI_sCell{1,1})),strcat(SaveMODDEMOFilePath,'Img_illus_lThreshPerim_PreROI_',timestamp,'.tif'));
        imwrite(uint8(255.*Img_Thresh_PreROI_sCell{2,1}),strcat(SaveMODDEMOFilePath,'Img_illus_fThresh_PreROI_',timestamp,'.tif'));
        imwrite(uint8(255.*bwperim(Img_Thresh_PreROI_sCell{2,1})),strcat(SaveMODDEMOFilePath,'Img_illus_fThreshPerim_PreROI_',timestamp,'.tif'));
        imwrite(uint8(255.*Img_Skel_PreROI_sCell{1,1}),strcat(SaveMODDEMOFilePath,'Img_illus_lSkel_PreROI_',timestamp,'.tif'));
        imwrite(uint8(255.*Img_Skel_PreROI_sCell{2,1}),strcat(SaveMODDEMOFilePath,'Img_illus_fSkel_PreROI_',timestamp,'.tif'));
        imwrite(uint8(255.*bwperim(Img_ROI)),strcat(SaveMODDEMOFilePath,'Img_illus_ROIPreROI_',timestamp,'.tif'));
    end

    clearvars Img_ROIDraw; 
    clearvars PolygonROI;
    clearvars FigHandle;
    clearvars SaveMODDEMOFilePath;


    %% Update ImgStack, Mask, Skeleton based on drawn ROI

    fprintf('Updating ImgStack, Mask, Skeleton based on drawn ROI...\n');

    [ImgStack_Crop,...
        Img_Thresh_ROI_Crop_sCell,rp_Thresh_ROI_Crop_sCell,...
        Img_Skel_ROI_Crop_sCell,rp_Skel_ROI_Crop_sCell,...
        Img_ROI_Crop,Img_Zstdev_Crop,...
        PolygonROI_startXY,ImgStack_ROICrop_cuboid,Img_ROI_Colocal_illus,...
        ImgStack_PreROICrop_Sz] = ...
        fPrepROIVesselImgStackParam(ImgStack,Img_ROI,Img_Zstdev,Img_Thresh_PreROI_sCell,Img_Skel_PreROI_sCell,PolygonROI_startXY,TimeParam);

    imwrite(Img_ROI_Colocal_illus,strcat(SaveMODFilePath,'ImgThresh_ImgROI_Colocalize.tif'));
    close all;

    clearvars Img_ROI_Colocal_illus;


    %% Collect to-be-saved Cropped ROI Information in structure

    oROI = struct;

    oROI.Img_Thresh_sCell = Img_Thresh_ROI_Crop_sCell; 
    oROI.rp_Thresh_sCell = rp_Thresh_ROI_Crop_sCell;

    oROI.Img_Skel_sCell = Img_Skel_ROI_Crop_sCell; 
    oROI.rp_Skel_sCell = rp_Skel_ROI_Crop_sCell;

    oROI.Img_Mask = Img_ROI_Crop; 
    oROI.Img_Zstdev = Img_Zstdev_Crop;

    oROI.ROICrop_cuboid = ImgStack_ROICrop_cuboid;

    oROI.PolygonROI_startXY = PolygonROI_startXY;

    oROI.ImgStack_PreROICrop_Sz = ImgStack_PreROICrop_Sz; 

    clearvars Img_Thresh_ROI_Crop_sCell rp_Thresh_ROI_Crop_sCell;
    clearvars Img_Skel_ROI_Crop_sCell rp_Skel_ROI_Crop_sCell;
    clearvars Img_ROI_Crop Img_Zstdev_Crop;
    clearvars ImgStack_ROICrop_cuboid;
    clearvars PolygonROI_startXY;
    clearvars ImgStack_PreROICrop_Sz;


    %% Clear variables likely no longer relevant after this module

    clearvars ImgStack Img_Zstdev Img_ROI; 
    clearvars Img_Thresh_PreROI_PreClean_sCell;
    clearvars Img_Thresh_PreROI_sCell;
    clearvars Img_Skel_PreROI_sCell; 


    %% Save Workspace

    if ~MOD00_AutoRun_Option
        fprintf('Saving Workspace...\n');
        varList_MOD00 = who;
        varList_All = varList_MOD00;
        save(strcat(SaveMODFilePath,'Workspace_STEP01Var_',SampleIDString,'_',timestamp,'.mat'));
    end

end


%% Save script in directory

ScriptName=mfilename;
PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveMODFilePath);
publish(strcat(ScriptName,'.m'),PublishOptions);


