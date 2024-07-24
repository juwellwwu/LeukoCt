clc; close all hidden; fclose('all');
% clearvars;


%% Collect to-be-saved Cropped ROI Information in structure, if not done in previous module

if  ~logical(exist('varList_MOD00')) && ~logical(exist('varList_All'))
    varList_MOD00 = who;
    varList_All = varList_MOD00;
end

if ~logical(exist('oROI'))

    oROI = struct;

    oROI.Img_Thresh_sCell = Img_Thresh_ROI_Crop_sCell; 
    oROI.rp_Thresh_sCell = rp_Thresh_ROI_Crop_sCell;

    oROI.Img_Skel_sCell = Img_Skel_ROI_Crop_sCell; 
    oROI.rp_Skel_sCell = rp_Skel_ROI_Crop_sCell;

    oROI.Img_Mask = Img_ROI_Crop; 
    oROI.Img_Zstdev = Img_Zstdev_Crop;

    oROI.ROICrop_cuboid = ImgStack_ROICrop_cuboid;

    oROI.PolygonROI_startXY = PolygonROI_startXY;

end

if ~isfield(TimeParam, 'tMidPt_FrLength')
    tFr_HistEdge = horzcat(TimeParam.tMidPt_FrStart,TimeParam.tMidPt_FrEnd(end));
    [TimeParam.tMidPt_FrLength,~,~] = histcounts(1:1:TimeParam.tMidPt_FrEnd(end),tFr_HistEdge);
end

if ~logical(exist('CBC_FilenameString')) 
    CBC_FilenameString = '/Users/j/Documents/MATLAB/Data_In/CBC_In.xlsx';
end

if ~exist('DEMOOption')
    DEMOOption = 'n';
end

clearvars Img_Thresh_ROI_Crop_sCell rp_Thresh_ROI_Crop_sCell;
clearvars Img_Skel_ROI_Crop_sCell rp_Skel_ROI_Crop_sCell;
clearvars Img_ROI_Crop Img_Zstdev_Crop;
clearvars ImgStack_ROICrop_cuboid;
clearvars PolygonROI_startXY;

clearvars Print_Seg_Preprocess_Montage; 
clearvars tFr_HistEdge;


%% Create Output Folder

SaveMODFilePath = strcat(SaveFilePath,'STEP02_SCS_VesselDiameter/');
if exist(SaveMODFilePath,'dir')==7
    rmdir(SaveMODFilePath,'s'); 
end
mkdir(SaveMODFilePath);

if DEMOOption == 'y'
    SaveMODDEMOFilePath = strcat(SaveMODFilePath,'DEMOIllus/');
    mkdir(SaveMODDEMOFilePath);
end


%% Set cell diameter unit length in pixels

CellDiameter_px = 5/XY_PxLength;


%% Specify Img_Thresh_Crop & Img_Skel_Crop to use

Img_Thresh_Crop = oROI.Img_Thresh_sCell{2,1}; 
rp_Thresh_Crop = oROI.rp_Thresh_sCell{2,1};

Img_ThreshVL_Crop = oROI.Img_Thresh_sCell{1,1}; 
rp_ThreshVL_Crop = oROI.rp_Thresh_sCell{1,1};

Img_Skel_Crop = oROI.Img_Skel_sCell{2,1}; 
rp_Skel_Crop = oROI.rp_Skel_sCell{2,1};

Img_SkelVL_Crop = oROI.Img_Skel_sCell{1,1};
rp_SkelVL_Crop = oROI.rp_Skel_sCell{1,1};

NumSeg = TimeParam.Num_tSeg;
NumMidPt = TimeParam.Num_tMidPt;


%% Prepare Img_Thresh for Vessel Flow Velocity Calculation

Img_Thresh_Crop_Cell = cell(1,NumSeg,1);
if Vessel_Diameter_TimeDependence == 'n'
    for SegCt = 1:NumSeg
        Img_Thresh_Crop_Cell{1,SegCt,1} = Img_Thresh_Crop;
    end
    NumThreshPerim = 3; 
else
    NumThreshPerim = 1+2*NumSeg;
end


%% Prepare Img_Thresh for Vessel Diameter Calculation

Img_ThreshVL_Crop_Cell = cell(1,NumSeg,1);
if Vessel_Diameter_TimeDependence == 'n'
    for SegCt = 1:NumSeg
        Img_ThreshVL_Crop_Cell{1,SegCt,1} = Img_ThreshVL_Crop;
    end
    NumThreshVLPerim = 3; 
else
    NumThreshVLPerim = 1+2*NumSeg;
end


%% Extract Img_Thresh_Crop_Cell Perimeters for Vessel Flow Velocity Measurement

% == Prepare cell to hold perimeters in Img_Thresh_Crop_Cell
ThreshPerim_Crop_xy_tpCell = cell(NumThreshPerim,1);
rp_ThreshPerim_Crop_tpCell = cell(NumThreshPerim,1); 
ImgStack_ThreshPerim_Crop = ... 
    zeros(size(Img_Thresh_Crop,1),size(Img_Thresh_Crop,2),1+NumSeg*2,'logical');

% == Load first element = Skeleton 
ImgStack_ThreshPerim_Crop(:,:,1) = Img_Skel_Crop; 
rp_ThreshPerim_Crop_tpCell{1,1} = regionprops(Img_Skel_Crop,'Area','PixelList','PixelIdxList');
ThreshPerim_Crop_xy_tpCell{1,1} = rp_ThreshPerim_Crop_tpCell{1,1}.PixelList;

% = Load other elements from Img_Thresh_Crop_Cell
for PerimPairCt = 1:(0.5*(NumThreshPerim-1))
        
        [Img_ThreshPerim_pCell,rp_ThreshPerim_pCell] = ... 
            fCreateThreshPerim(Img_Thresh_Crop_Cell{1,PerimPairCt,1},oROI.Img_Mask);

        % Load 1st perimeter (Vessel Boundary 1)
        ImgStack_ThreshPerim_Crop(:,:,PerimPairCt*2) = Img_ThreshPerim_pCell{1,1};
        rp_ThreshPerim_Crop_tpCell{PerimPairCt*2,1} = rp_ThreshPerim_pCell{1,1};
        ThreshPerim_Crop_xy_tpCell{PerimPairCt*2,1} = rp_ThreshPerim_Crop_tpCell{PerimPairCt*2,1}.PixelList;

        % Load 2nd perimeter (Vessel Boundary 2)
        ImgStack_ThreshPerim_Crop(:,:,PerimPairCt*2+1) = Img_ThreshPerim_pCell{2,1};
        rp_ThreshPerim_Crop_tpCell{PerimPairCt*2+1,1} = rp_ThreshPerim_pCell{2,1};
        ThreshPerim_Crop_xy_tpCell{PerimPairCt*2+1,1} = rp_ThreshPerim_Crop_tpCell{PerimPairCt*2+1,1}.PixelList;

end

% % % figure; 
% % % imshow(max(ImgStack_ThreshPerim_Crop,[],3));

clearvars Img_ThreshPerim_pCell rp_ThreshPerim_pCell PerimPairCt;
clearvars PerimPairCt;


%% Extract Img_Thresh_Crop_Cell Perimeters for Vessel Diameter Measurement

% == Prepare cell to hold perimeters in Img_Thresh_Crop_Cell
ThreshVLPerim_Crop_xy_tpCell = cell(3,1);
rp_ThreshVLPerim_Crop_tpCell = cell(3,1); 
ImgStack_ThreshVLPerim_Crop = ... 
    zeros(size(Img_ThreshVL_Crop,1),size(Img_ThreshVL_Crop,2),1+NumSeg*2,'logical');

% == Load first element = Skeleton 
ImgStack_ThreshVLPerim_Crop(:,:,1) = Img_SkelVL_Crop; 
rp_ThreshVLPerim_Crop_tpCell{1,1} = regionprops(Img_SkelVL_Crop,'Area','PixelList','PixelIdxList');
ThreshVLPerim_Crop_xy_tpCell{1,1} = rp_ThreshVLPerim_Crop_tpCell{1,1}.PixelList;

% = Load other elements from Img_Thresh_Crop_Cell
for PerimPairCt = 1:(0.5*(NumThreshVLPerim-1))
        
        [Img_ThreshVLPerim_pCell,rp_ThreshVLPerim_pCell] = ...
            fCreateThreshPerim(Img_ThreshVL_Crop_Cell{1,PerimPairCt,1},oROI.Img_Mask);

        % Load 1st perimeter (Vessel Boundary 1)
        ImgStack_ThreshVLPerim_Crop(:,:,PerimPairCt*2) = Img_ThreshVLPerim_pCell{1,1};
        rp_ThreshVLPerim_Crop_tpCell{PerimPairCt*2,1} = rp_ThreshVLPerim_pCell{1,1};
        ThreshVLPerim_Crop_xy_tpCell{PerimPairCt*2,1} = rp_ThreshVLPerim_Crop_tpCell{PerimPairCt*2,1}.PixelList;

        % Load 2nd perimeter (Vessel Boundary 2)
        ImgStack_ThreshVLPerim_Crop(:,:,PerimPairCt*2+1) = Img_ThreshVLPerim_pCell{2,1};
        rp_ThreshVLPerim_Crop_tpCell{PerimPairCt*2+1,1} = rp_ThreshVLPerim_pCell{2,1};
        ThreshVLPerim_Crop_xy_tpCell{PerimPairCt*2+1,1} = rp_ThreshVLPerim_Crop_tpCell{PerimPairCt*2+1,1}.PixelList;

end

% % % figure; 
% % % imshow(max(ImgStack_ThreshPerim_Crop,[],3));

clearvars Img_ThreshVLPerim_pCell rp_ThreshVLPerim_pCell PerimPairCt;
clearvars PerimPairCt;


%% Order pixels of Img_Thresh Perimeters from one endpt to the other

for ThreshPerimCt = 1:NumThreshPerim 

    ThreshPerimPxList = ThreshPerim_Crop_xy_tpCell{ThreshPerimCt,1};

    Img_Endpt_Temp = bwmorph(ImgStack_ThreshPerim_Crop(:,:,ThreshPerimCt),'endpoints');
    rp_Endpt_Temp = regionprops(Img_Endpt_Temp,'PixelList','PixelIdxList');

    % Determine the endpoint closer to PolygonROI_startXY
    [~, PolygonROI_startXY_EndptIdx]=pdist2(cat(1,rp_Endpt_Temp.PixelList),oROI.PolygonROI_startXY,'euclidean','Smallest',1);
    [~, Endpt_Loc] = ismember(rp_Endpt_Temp(PolygonROI_startXY_EndptIdx).PixelList,ThreshPerimPxList,'rows');

    % Order Skeleton pixels
    ThreshPerim_Crop_xy_tpCell{ThreshPerimCt,1} = fOrderSkelPx(ThreshPerimPxList,Endpt_Loc);

end

clearvars ThreshPerimCt ThreshPerimPxList Img_Endpt_Temp rp_Endpt_Temp;
clearvars PolygonROI_startXY_EndptIdx Endpt_Loc;


%% Order pixels of Img_Thresh Perimeters (for Vessel Diameter Measurement) from one endpt to the other

for ThreshPerimCt = 1:NumThreshVLPerim 

    ThreshVLPerimPxList = ThreshVLPerim_Crop_xy_tpCell{ThreshPerimCt,1};

    Img_Endpt_Temp = bwmorph(ImgStack_ThreshVLPerim_Crop(:,:,ThreshPerimCt),'endpoints');
    rp_Endpt_Temp = regionprops(Img_Endpt_Temp,'PixelList','PixelIdxList');

    % Determine the endpoint closer to PolygonROI_startXY
    [~, PolygonROI_startXY_EndptIdx]=pdist2(cat(1,rp_Endpt_Temp.PixelList),oROI.PolygonROI_startXY,'euclidean','Smallest',1);
    [~, Endpt_Loc] = ismember(rp_Endpt_Temp(PolygonROI_startXY_EndptIdx).PixelList,ThreshVLPerimPxList,'rows');

    % Order Skeleton pixels
    ThreshVLPerim_Crop_xy_tpCell{ThreshPerimCt,1} = fOrderSkelPx(ThreshVLPerimPxList,Endpt_Loc);

end

clearvars ThreshPerimCt ThreshVLPerimPxList Img_Endpt_Temp rp_Endpt_Temp;
clearvars PolygonROI_startXY_EndptIdx Endpt_Loc;


%% Find Normal Lines to Original Skeleton (Skeleton defined for Vessel Flow Velocity measurement)

% == Find Normal Lines to Original Skeleton
[SkelNorm_xlim,SkelNorm_ylim,SkelNorm_FullLn_Endpt_RowIdx,Skel_Crop_Curvature_px] = ...
    fCreateSkelNormal(ThreshPerim_Crop_xy_tpCell{1,1},ceil(CellDiameter_px*1.00),0.5*(max(size(Img_Skel_Crop))));

% == Update original (ordered) skeleton
ThreshPerim_Crop_xy_tpCell{1,1}((SkelNorm_FullLn_Endpt_RowIdx(2)+1):end,:) = []; 
ThreshPerim_Crop_xy_tpCell{1,1}(1:(SkelNorm_FullLn_Endpt_RowIdx(1)-1),:) = [];

Img_Temp = zeros(size(Img_Skel_Crop),'logical');
Img_Temp(sub2ind(size(Img_Skel_Crop),ThreshPerim_Crop_xy_tpCell{1,1}(:,2),ThreshPerim_Crop_xy_tpCell{1,1}(:,1))) = 1;
ImgStack_ThreshPerim_Crop(:,:,1) = Img_Temp;

rp_ThreshPerim_Crop_tpCell{1,1} = regionprops(ImgStack_ThreshPerim_Crop(:,:,1),'Area','PixelList','PixelIdxList');

% Img_Skel_Crop = ImgStack_ThreshPerim_Crop(:,:,1);
% rp_Skel_Crop = rp_ThreshPerim_Crop_tpCell{1,1}; 

clearvars Img_Temp SkelNorm_FullLn_Endpt_RowIdx;


%% Find Normal Lines to Original Skeleton (Skeleton defined for Vessel Diameter measurement)

% == Find Normal Lines to Original Skeleton
[SkelVLNorm_xlim,SkelVLNorm_ylim,SkelVLNorm_FullLn_Endpt_RowIdx,~] = ...
    fCreateSkelNormal(ThreshVLPerim_Crop_xy_tpCell{1,1},ceil(CellDiameter_px*1.00),0.5*(max(size(Img_SkelVL_Crop))));

% == Update original (ordered) skeleton
ThreshVLPerim_Crop_xy_tpCell{1,1}((SkelVLNorm_FullLn_Endpt_RowIdx(2)+1):end,:) = []; % Must clear this end first
ThreshVLPerim_Crop_xy_tpCell{1,1}(1:(SkelVLNorm_FullLn_Endpt_RowIdx(1)-1),:) = [];

Img_Temp = zeros(size(Img_SkelVL_Crop),'logical');
Img_Temp(sub2ind(size(Img_SkelVL_Crop),ThreshVLPerim_Crop_xy_tpCell{1,1}(:,2),ThreshVLPerim_Crop_xy_tpCell{1,1}(:,1))) = 1;
ImgStack_ThreshVLPerim_Crop(:,:,1) = Img_Temp;

rp_ThreshVLPerim_Crop_tpCell{1,1} = regionprops(ImgStack_ThreshVLPerim_Crop(:,:,1),'Area','PixelList','PixelIdxList');

% Img_SkelVL_Crop = ImgStack_ThreshVLPerim_Crop(:,:,1);
% rp_SkelVL_Crop = rp_ThreshVLPerim_Crop_tpCell{1,1}; 

clearvars Img_Temp SkelVLNorm_FullLn_Endpt_RowIdx;


%% Calculating Approx Vessel Diameter, Pre Vessel Block Definition

Img_DistMap_Crop_tpCell = cell(1,NumThreshVLPerim,1);
Vessel_Diameter_PreBlkDefine_px = zeros(height(ThreshVLPerim_Crop_xy_tpCell{1,1}),NumThreshVLPerim);

for tpCt = 1:NumThreshVLPerim
    [Img_DistMap_Crop_tpCell{1,tpCt,1},~] = ...
        bwdist(bwperim(ImgStack_ThreshVLPerim_Crop(:,:,tpCt)),'euclidean');
    Vessel_Diameter_PreBlkDefine_px(:,tpCt) = ...
        Img_DistMap_Crop_tpCell{1,tpCt,1}(sub2ind(size(Img_Skel_Crop),ThreshVLPerim_Crop_xy_tpCell{1,1}(:,2),ThreshVLPerim_Crop_xy_tpCell{1,1}(:,1)));
end

Vessel_Diameter_PreBlkDefine_px = Vessel_Diameter_PreBlkDefine_px(:,2:2:end-1) + Vessel_Diameter_PreBlkDefine_px(:,3:2:end);
Vessel_Diameter_Max_PreBlkDefine_px = max(Vessel_Diameter_PreBlkDefine_px,[],2); 

[Vessel_Diameter_PreBlkDefine_px_MapPxDist,Vessel_Diameter_PreBlkDefine_px_MapIdx] = ...
    pdist2(ThreshVLPerim_Crop_xy_tpCell{1,1},ThreshPerim_Crop_xy_tpCell{1,1},'euclidean','Smallest',1);

clearvars Img_DistMap_Crop_tpCell;
clearvars tpCt;


%% Apply Normal Lines to Original Skeleton (Skeleton defined for Vessel Flow Velocity measurement) to parameters for Vessel Diameter measurement 

SkelVLNorm_xlim = SkelNorm_xlim;
SkelVLNorm_ylim = SkelNorm_ylim;

ThreshVLPerim_Crop_xy_tpCell{1,1} = ThreshPerim_Crop_xy_tpCell{1,1};
ImgStack_ThreshVLPerim_Crop(:,:,1) = ImgStack_ThreshPerim_Crop(:,:,1);
rp_ThreshVLPerim_Crop_tpCell{1,1} = rp_ThreshPerim_Crop_tpCell{1,1}; 


%% Determine intersections of original skeleton normals w/ all skeleton dilations + Img_Thresh perimeters, for Vessel Flow Velocity measurement 

% Prepare cells to hold Intersection information for all threshold
% perimeters
Intersect_ThreshPerim_xy_tpCell = cell(NumThreshPerim,height(ThreshPerim_Crop_xy_tpCell{1,1}));
Intersect_ThreshPerimPxIdx_tpCell = cell(NumThreshPerim,height(ThreshPerim_Crop_xy_tpCell{1,1}));

for ThreshPerimPxCt = 1:height(ThreshPerim_Crop_xy_tpCell{1,1})

    % = Create image of line normal of original skeleton @ pixel SkelPxCt
    Img_Norm_Temp = fCreateLineROIImg(SkelNorm_xlim(ThreshPerimPxCt,:),SkelNorm_ylim(ThreshPerimPxCt,:),width(ImgStack_Crop),height(ImgStack_Crop));

    % = Create and apply mask    
    Img_Norm_Mask_Temp = zeros(size(Img_Norm_Temp),'logical');
    Img_Norm_Mask_Temp(sub2ind(size(Img_Norm_Temp),ThreshPerim_Crop_xy_tpCell{1,1}(ThreshPerimPxCt,2),ThreshPerim_Crop_xy_tpCell{1,1}(ThreshPerimPxCt,1))) = 1;
    Img_Norm_Mask_Temp = ... 
        imdilate(Img_Norm_Mask_Temp,...
        strel('disk',double(ceil(3+Vessel_Diameter_PreBlkDefine_px_MapPxDist(ThreshPerimPxCt)+0.5*Vessel_Diameter_Max_PreBlkDefine_px(Vessel_Diameter_PreBlkDefine_px_MapIdx(ThreshPerimPxCt))))));

    Img_Norm_Temp = Img_Norm_Temp.*Img_Norm_Mask_Temp;

    % = Find intersects to all threshold perimeters by skeleton normal Line
    [Intersect_ThreshPerimPxIdx_tpCell,Intersect_ThreshPerim_xy_tpCell,~,~] = ...
        fFindSkelNormIntersect(Img_Norm_Temp,ThreshPerim_Crop_xy_tpCell,ThreshPerimPxCt,Intersect_ThreshPerimPxIdx_tpCell,Intersect_ThreshPerim_xy_tpCell);
    
end

close all;

clearvars SkelNorm_xlim SkelNorm_ylim;
clearvars Img_Norm_Temp Img_Norm_Mask_Temp
clearvars ThreshPerimPxCt;


%% Determine intersections of original skeleton normals w/ all skeleton dilations + Img_Thresh perimeters, for Vessel Diameter measurement 

% Prepare cells to hold Intersection information for all threshold
% perimeters
Intersect_ThreshVLPerim_xy_tpCell = cell(NumThreshVLPerim,height(ThreshVLPerim_Crop_xy_tpCell{1,1}));
Intersect_ThreshVLPerimPxIdx_tpCell = cell(NumThreshVLPerim,height(ThreshVLPerim_Crop_xy_tpCell{1,1}));

for ThreshPerimPxCt = 1:height(ThreshVLPerim_Crop_xy_tpCell{1,1})

    % = Create image of line normal of original skeleton @ pixel SkelPxCt
    Img_Norm_Temp = fCreateLineROIImg(SkelVLNorm_xlim(ThreshPerimPxCt,:),SkelVLNorm_ylim(ThreshPerimPxCt,:),width(ImgStack_Crop),height(ImgStack_Crop));

    % = Create and apply mask   
    Img_Norm_Mask_Temp = zeros(size(Img_Norm_Temp),'logical');
    Img_Norm_Mask_Temp(sub2ind(size(Img_Norm_Temp),ThreshVLPerim_Crop_xy_tpCell{1,1}(ThreshPerimPxCt,2),ThreshVLPerim_Crop_xy_tpCell{1,1}(ThreshPerimPxCt,1))) = 1 ;
    Img_Norm_Mask_Temp = ...
        imdilate(Img_Norm_Mask_Temp,...
        strel('disk',double(ceil(3+Vessel_Diameter_PreBlkDefine_px_MapPxDist(ThreshPerimPxCt)+0.5*Vessel_Diameter_Max_PreBlkDefine_px(Vessel_Diameter_PreBlkDefine_px_MapIdx(ThreshPerimPxCt))))));

    Img_Norm_Temp = Img_Norm_Temp.*Img_Norm_Mask_Temp;
    
    % = Find intersects to all threshold perimeters by skeleton normal Line
    [Intersect_ThreshVLPerimPxIdx_tpCell,Intersect_ThreshVLPerim_xy_tpCell,~,~] = ...
        fFindSkelNormIntersect(Img_Norm_Temp,ThreshVLPerim_Crop_xy_tpCell,ThreshPerimPxCt,Intersect_ThreshVLPerimPxIdx_tpCell,Intersect_ThreshVLPerim_xy_tpCell);

end

close all;

clearvars SkelVLNorm_xlim SkelVLNorm_ylim;
clearvars Img_Norm_Temp Img_Norm_Mask_Temp
clearvars ThreshPerimPxCt;


%% Use intersects to equalize length of all time-dependent Img_Thresh perimeters (for Vessel Flow Velocity measurement)

[ThreshPerim_Crop_xy_tpCell,ThreshPerim_Crop_xy_Flag_tpCell,ImgStack_ThreshPerim_Crop,rp_ThreshPerim_Crop_tpCell] = ...
    fMatchSkelLength(ThreshPerim_Crop_xy_tpCell,Intersect_ThreshPerimPxIdx_tpCell,ImgStack_ThreshPerim_Crop,rp_ThreshPerim_Crop_tpCell);

clearvars Intersect_ThreshPerimPxIdx_tpCell Intersect_ThreshPerim_xy_tpCell;


%% Use intersects to equalize length of all time-dependent Img_Thresh perimeters (for Vessel Diameter measurement)

[ThreshVLPerim_Crop_xy_tpCell,ThreshVLPerim_Crop_xy_Flag_tpCell,ImgStack_ThreshVLPerim_Crop,rp_ThreshVLPerim_Crop_tpCell] = ...
    fMatchSkelLength(ThreshVLPerim_Crop_xy_tpCell,Intersect_ThreshVLPerimPxIdx_tpCell,ImgStack_ThreshVLPerim_Crop,rp_ThreshVLPerim_Crop_tpCell);

clearvars Intersect_ThreshVLPerimPxIdx_tpCell Intersect_ThreshVLPerim_xy_tpCell;


%% Finalize ThreshPerim: Divide Img_Thresh Perimeters into Vessel Blocks (for Vessel Flow Velocity measurement)

% = Determine length of shortest Img_Thresh perimeter; will decide # vessel blocks
MinThreshPerimLength = min(min(cell2mat(cellfun(@height,ThreshPerim_Crop_xy_tpCell(1:NumThreshPerim),'UniformOutput',false)),[],1), ...
    min(cell2mat(cellfun(@height,ThreshVLPerim_Crop_xy_tpCell(1:NumThreshPerim),'UniformOutput',false)),[],1)); 

% = Divide Img_Thresh perimeters into blocks of x cell-diameter length blocks
SkelBlockLength = unique(round(CellDiameter_px*(max(5,(NumCellDiameter_in_SkelBlockLength-2)):0.1:(NumCellDiameter_in_SkelBlockLength+2))));
SkelBlockLength_Rem = rem(MinThreshPerimLength,SkelBlockLength); 
[~,SkelBlockLength_Rem_MinIdx] = min(SkelBlockLength_Rem); 
SkelBlockLength =  SkelBlockLength(SkelBlockLength_Rem_MinIdx);
NumBlk = floor(MinThreshPerimLength/SkelBlockLength);

% = Divide Img_Thresh perimenters into vessel blocks
[ThreshPerim_Crop_xy_tpCell,ThreshPerim_Crop_xy_Flag_tpCell] = ...
    fDefineVesselBlock(NumBlk,SkelBlockLength,ThreshPerim_Crop_xy_tpCell,ThreshPerim_Crop_xy_Flag_tpCell);

clearvars MinThreshPerimLength SkelBlockLength_Rem SkelBlockLength_Rem_MinIdx;
clearvars ImgStack_ThreshPerim_Crop rp_ThreshPerim_Crop_Cell; 
clearvars *PreBlkDefine*; 
clearvars Skel_Crop_Curvature_px; 


%% Finalize ThreshPerim: Divide Img_Thresh Perimeters into Vessel Blocks (for Vessel Diameter measurement)

[ThreshVLPerim_Crop_xy_tpCell,ThreshVLPerim_Crop_xy_Flag_tpCell] = ...
    fDefineVesselBlock(NumBlk,SkelBlockLength,ThreshVLPerim_Crop_xy_tpCell,ThreshVLPerim_Crop_xy_Flag_tpCell);

clearvars ImgStack_ThreshVLPerim_Crop rp_ThreshVLPerim_Crop_Cell; 
clearvars *PreBlkDefine*; 


%%  Create Vessel Flow Profile Radius Cell, for SCS creation

[fRadiusL,fDiameter] = fCalcVesselSize(ThreshPerim_Crop_xy_tpCell); 

[tMidPt_SkPx_fRadiusL,~] = fConvertMtx_tSegtoMidPt(fRadiusL,TimeParam); 
[tMidPt_SkPx_fDiameter,~] = fConvertMtx_tSegtoMidPt(fDiameter,TimeParam); 

tMidPtwMean_SkPx_fDiameter = ... 
    sum(repmat(TimeParam.tMidPt_FrLength,size(tMidPt_SkPx_fDiameter,1),1,1).*tMidPt_SkPx_fDiameter,2)./sum(repmat(TimeParam.tMidPt_FrLength,size(tMidPt_SkPx_fDiameter,1),1,1),2);
tMidPtMin_SkPx_fDiameter = min(tMidPt_SkPx_fDiameter,[],2); 
tMidPtMax_SkPx_fDiameter = max(tMidPt_SkPx_fDiameter,[],2); 


%%  Create Vessel Lumen Radius Cell, for Vessel Diameter Measurement

[lRadiusL,lDiameter] = fCalcVesselSize(ThreshVLPerim_Crop_xy_tpCell); 

[tMidPt_SkPx_lRadiusL] = fConvertMtx_tSegtoMidPt(lRadiusL,TimeParam); 
[tMidPt_SkPx_lDiameter, ~] = fConvertMtx_tSegtoMidPt(lDiameter,TimeParam); 

tMidPtwMean_SkPx_lDiameter = ... 
    sum(repmat(TimeParam.tMidPt_FrLength,size(tMidPt_SkPx_lDiameter,1),1,1).*tMidPt_SkPx_lDiameter,2)./sum(repmat(TimeParam.tMidPt_FrLength,size(tMidPt_SkPx_lDiameter,1),1,1),2);
tMidPtMin_SkPx_lDiameter = min(tMidPt_SkPx_lDiameter,[],2); 
tMidPtMax_SkPx_lDiameter = max(tMidPt_SkPx_lDiameter,[],2); 


%%  Plot Vessel Lumen Diameter

fprintf('Plotting Vessel Lumen Diameter...\n');

FigHandle0 = figure('Position',[10 500 1800 300]); % [left bottom width height]
hold on;
errorbar(1:height(tMidPtwMean_SkPx_lDiameter),tMidPtwMean_SkPx_lDiameter.*XY_PxLength,...
    (tMidPtwMean_SkPx_lDiameter-tMidPtMin_SkPx_lDiameter).*XY_PxLength,... 
    (tMidPtMax_SkPx_lDiameter-tMidPtwMean_SkPx_lDiameter).*XY_PxLength,...
    'Color',[1 0 0], 'LineWidth',1,"CapSize",3);
xlim([1 height(tMidPtwMean_SkPx_lDiameter)]);
xticks(0:5:5*ceil(height(tMidPtwMean_SkPx_lDiameter)/5));
xlabel('Pixel Along Skeleton');
ylabel('Vessel Diameter (um)');
title({strcat('Vessel Lumen Diameter (um): Mean, Min (Lower Error Bar), Max (Upper Error Bar) of all tMidPts')...
    strcat('Overall Mean',32,num2str(mean(tMidPtwMean_SkPx_lDiameter,'all').*XY_PxLength,'%0.2f'),';',32,...
    'Min:',32,num2str(min(tMidPtMin_SkPx_lDiameter,[],'all').*XY_PxLength,'%0.2f'),'; Max:',32,num2str(max(tMidPtMax_SkPx_lDiameter,[],'all').*XY_PxLength,'%0.2f')),...
    strcat('By Vessel Blk Mean (um):',32,char(strjoin(string(splitapply(@mean,tMidPtwMean_SkPx_lDiameter,(repelem((1:1:NumBlk),SkelBlockLength))')))),';',32,...
    'Min:',32,char(strjoin(string(splitapply(@min,tMidPtMin_SkPx_lDiameter,(repelem((1:1:NumBlk),SkelBlockLength))')))),...
    '; Max:',32,char(strjoin(string(splitapply(@max,tMidPtMax_SkPx_lDiameter,(repelem((1:1:NumBlk),SkelBlockLength))'))))),...
    strcat('SkelBlockLength (px):',32,num2str(SkelBlockLength),'; (um):',32,num2str(SkelBlockLength.*XY_PxLength))},...
    'FontSize',8,'FontWeight','Normal','interpreter','none');

% = Save Plot
exportgraphics(FigHandle0,strcat(SaveMODFilePath,'tMidPt_Vessel_Lumen_Diameter_um_Plot_',timestamp,'.tif'),'Resolution',150);
% close all;

fprintf(strcat('\n... Mean Vessel ROI Lumen Diameter (um) =',32,num2str(mean(tMidPt_SkPx_lDiameter.*XY_PxLength,'all'),'%0.2f')));
fprintf(strcat('\n... Min Vessel ROI Lumen Diameter (um) =',32,num2str(min(max(tMidPt_SkPx_lDiameter,[],2),[],'all').*XY_PxLength,'%0.2f')));
fprintf(strcat('\n... Max Vessel ROI Lumen Diameter (um) =',32,num2str(max(max(tMidPt_SkPx_lDiameter,[],2),[],'all').*XY_PxLength,'%0.2f')));
fprintf('\n\n');

clearvars Vessel_Diameter_PreBlkDefine_px  Vessel_Diameter_Max_PreBlkDefine_px; % Obsolete
clearvars FigHandle0*;


%% Calculate NumSkel_VsWidth and NumSkel from Vessel Size

tMidPt_fRadius_Max = max(max(tMidPt_SkPx_fRadiusL,[],2),[],1); 
tMidPt_fRadius_Min = min(max(tMidPt_SkPx_fRadiusL,[],2),[],1); 

NumSkel_VsWidth = 1+2*ceil(tMidPt_fRadius_Max/(1.0/XY_PxLength)); 
NumSkel = max((1+2*ceil(0.20*tMidPt_fRadius_Min/(1.0/XY_PxLength))),3); 

clearvars tMidPt_fRadius_Max tMidPt_fRadius_Min;


%% Calculate Vessel Flow Profile Area (px^2)

% = Calculate as matrix first
tMidPt_SkPx_fDiameter_AllSk_Cell = tMidPt_SkPx_fDiameter; 
tMidPt_SkPx_fDiameter_AllSk_Cell = cell2mat(permute(mat2cell(tMidPt_SkPx_fDiameter_AllSk_Cell,repelem(SkelBlockLength,NumBlk),NumMidPt),[3 2 1])); 

tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell = ...  
    sum(repmat(TimeParam.tMidPt_FrLength,SkelBlockLength,1,NumBlk).*tMidPt_SkPx_fDiameter_AllSk_Cell,2)./sum(repmat(TimeParam.tMidPt_FrLength,SkelBlockLength,1,NumBlk),2);
tMidPtMax_SkPxMean_fDiameter_AllSk_Cell = max(tMidPt_SkPx_fDiameter_AllSk_Cell,[],2);

tMidPt_SkPx_fDiameter_AllSk_Cell = mean(tMidPt_SkPx_fDiameter_AllSk_Cell,1); 

tMidPtwMean_SkPxMin_fDiameter_AllSk_Cell = min(tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell,[],1); 
tMidPtwMean_SkPxMax_fDiameter_AllSk_Cell = max(tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell,[],1); 
tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell = mean(tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell,1); 

tMidPtMax_SkPxMin_fDiameter_AllSk_Cell = min(tMidPtMax_SkPxMean_fDiameter_AllSk_Cell,[],1); 
tMidPtMax_SkPxMax_fDiameter_AllSk_Cell = max(tMidPtMax_SkPxMean_fDiameter_AllSk_Cell,[],1); 
tMidPtMax_SkPxMean_fDiameter_AllSk_Cell = mean(tMidPtMax_SkPxMean_fDiameter_AllSk_Cell,1); 

tMidPt_fArea_AllSk_Cell = pi*(0.5*tMidPt_SkPx_fDiameter_AllSk_Cell).^2; 

tMidPtwMean_SkPxMean_fArea_AllSk_Cell = pi*(0.5*tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell).^2;  
tMidPtMax_SkPxMean_fArea_AllSk_Cell = pi*(0.5*tMidPtMax_SkPxMean_fDiameter_AllSk_Cell).^2;  

% = Convert to cell
tMidPt_fDiameter_Cell = repmat(tMidPt_SkPx_fDiameter_AllSk_Cell,NumSkel_VsWidth,1,1); 
tMidPt_fArea_Cell = repmat(tMidPt_fArea_AllSk_Cell,NumSkel_VsWidth,1,1); 

tMidPt_SkPx_fDiameter_AllSk_Cell = mat2cell(tMidPt_SkPx_fDiameter_AllSk_Cell,1,repelem(1,NumMidPt),repelem(1,NumBlk)); 

tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell = mat2cell(tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtwMean_SkPxMin_fDiameter_AllSk_Cell = mat2cell(tMidPtwMean_SkPxMin_fDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtwMean_SkPxMax_fDiameter_AllSk_Cell = mat2cell(tMidPtwMean_SkPxMax_fDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 

tMidPtMax_SkPxMean_fDiameter_AllSk_Cell = mat2cell(tMidPtMax_SkPxMean_fDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtMax_SkPxMin_fDiameter_AllSk_Cell = mat2cell(tMidPtMax_SkPxMin_fDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtMax_SkPxMax_fDiameter_AllSk_Cell = mat2cell(tMidPtMax_SkPxMax_fDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 

tMidPt_fArea_AllSk_Cell = mat2cell(tMidPt_fArea_AllSk_Cell,1,repelem(1,NumMidPt),repelem(1,NumBlk)); 

tMidPtwMean_SkPxMean_fArea_AllSk_Cell = mat2cell(tMidPtwMean_SkPxMean_fArea_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtMax_SkPxMean_fArea_AllSk_Cell = mat2cell(tMidPtMax_SkPxMean_fArea_AllSk_Cell,1,1,repelem(1,NumBlk)); 


%% Calculate Vessel Lumen Area (px^2)

% = Calculate as matrix first
tMidPt_SkPx_lDiameter_AllSk_Cell = tMidPt_SkPx_lDiameter;
tMidPt_SkPx_lDiameter_AllSk_Cell = cell2mat(permute(mat2cell(tMidPt_SkPx_lDiameter_AllSk_Cell,repelem(SkelBlockLength,NumBlk),NumMidPt),[3 2 1])); 

tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell = ...  
    sum(repmat(TimeParam.tMidPt_FrLength,SkelBlockLength,1,NumBlk).*tMidPt_SkPx_lDiameter_AllSk_Cell,2)./sum(repmat(TimeParam.tMidPt_FrLength,SkelBlockLength,1,NumBlk),2);
tMidPtMax_SkPxMean_lDiameter_AllSk_Cell = max(tMidPt_SkPx_lDiameter_AllSk_Cell,[],2);

tMidPt_SkPx_lDiameter_AllSk_Cell = mean(tMidPt_SkPx_lDiameter_AllSk_Cell,1);  

tMidPtwMean_SkPxMin_lDiameter_AllSk_Cell = min(tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell,[],1); 
tMidPtwMean_SkPxMax_lDiameter_AllSk_Cell = max(tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell,[],1); 
tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell = mean(tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell,1); 

tMidPtMax_SkPxMin_lDiameter_AllSk_Cell = min(tMidPtMax_SkPxMean_lDiameter_AllSk_Cell,[],1); 
tMidPtMax_SkPxMax_lDiameter_AllSk_Cell = max(tMidPtMax_SkPxMean_lDiameter_AllSk_Cell,[],1); 
tMidPtMax_SkPxMean_lDiameter_AllSk_Cell = mean(tMidPtMax_SkPxMean_lDiameter_AllSk_Cell,1); 

tMidPt_lArea_AllSk_Cell = pi*(0.5*tMidPt_SkPx_lDiameter_AllSk_Cell).^2; 

tMidPtwMean_SkPxMean_lArea_AllSk_Cell = pi*(0.5*tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell).^2; 
tMidPtMax_SkPxMean_lArea_AllSk_Cell = pi*(0.5*tMidPtMax_SkPxMean_lDiameter_AllSk_Cell).^2;

% = Convert to cell
tMidPt_lDiameter_Cell = repmat(tMidPt_SkPx_lDiameter_AllSk_Cell,NumSkel_VsWidth,1,1); 
tMidPt_lArea_Cell = repmat(tMidPt_lArea_AllSk_Cell,NumSkel_VsWidth,1,1); 

tMidPt_SkPx_lDiameter_AllSk_Cell = mat2cell(tMidPt_SkPx_lDiameter_AllSk_Cell,1,repelem(1,NumMidPt),repelem(1,NumBlk)); 

tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell = mat2cell(tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtwMean_SkPxMin_lDiameter_AllSk_Cell = mat2cell(tMidPtwMean_SkPxMin_lDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtwMean_SkPxMax_lDiameter_AllSk_Cell = mat2cell(tMidPtwMean_SkPxMax_lDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 

tMidPtMax_SkPxMean_lDiameter_AllSk_Cell = mat2cell(tMidPtMax_SkPxMean_lDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtMax_SkPxMin_lDiameter_AllSk_Cell = mat2cell(tMidPtMax_SkPxMin_lDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtMax_SkPxMax_lDiameter_AllSk_Cell = mat2cell(tMidPtMax_SkPxMax_lDiameter_AllSk_Cell,1,1,repelem(1,NumBlk)); 

tMidPt_lArea_AllSk_Cell = mat2cell(tMidPt_lArea_AllSk_Cell,1,repelem(1,NumMidPt),repelem(1,NumBlk)); 

tMidPtwMean_SkPxMean_lArea_AllSk_Cell = mat2cell(tMidPtwMean_SkPxMean_lArea_AllSk_Cell,1,1,repelem(1,NumBlk)); 
tMidPtMax_SkPxMean_lArea_AllSk_Cell = mat2cell(tMidPtMax_SkPxMean_lArea_AllSk_Cell,1,1,repelem(1,NumBlk)); 


%% Collect to-be-saved Vessel Radius and Diameter in structure

oRadius = struct;

% = Vessel Flow Profile Radius (px), Diameter (px), Area (px2)
oRadius.F.tMidPt_SkPx_RadiusL = tMidPt_SkPx_fRadiusL; 
oRadius.F.tMidPt_SkPx_Diameter = tMidPt_SkPx_fDiameter;
oRadius.F.tMidPtwMean_SkPx_Diameter = tMidPtwMean_SkPx_fDiameter;
oRadius.F.tMidPtMin_SkPx_Diameter = tMidPtMin_SkPx_fDiameter;
oRadius.F.tMidPtMax_SkPx_Diameter = tMidPtMax_SkPx_fDiameter;

oRadius.F.tMidPt_Diameter_Cell = tMidPt_fDiameter_Cell;
oRadius.F.tMidPt_AllSk_Diameter_Cell = tMidPt_SkPx_fDiameter_AllSk_Cell;
oRadius.F.tMidPtwMean_SkPxMean_AllSk_Diameter_Cell = tMidPtwMean_SkPxMean_fDiameter_AllSk_Cell;
oRadius.F.tMidPtwMean_SkPxMin_AllSk_Diameter_Cell = tMidPtwMean_SkPxMin_fDiameter_AllSk_Cell;
oRadius.F.tMidPtwMean_SkPxMax_AllSk_Diameter_Cell = tMidPtwMean_SkPxMax_fDiameter_AllSk_Cell;
oRadius.F.tMidPtMax_SkPxMean_AllSk_Diameter_Cell = tMidPtMax_SkPxMean_fDiameter_AllSk_Cell;
oRadius.F.tMidPtMax_SkPxMin_AllSk_Diameter_Cell = tMidPtMax_SkPxMin_fDiameter_AllSk_Cell;
oRadius.F.tMidPtMax_SkPxMax_AllSk_Diameter_Cell = tMidPtMax_SkPxMax_fDiameter_AllSk_Cell;

oRadius.F.tMidPt_Area_Cell = tMidPt_fArea_Cell;
oRadius.F.tMidPt_AllSk_Area_Cell = tMidPt_fArea_AllSk_Cell;
oRadius.F.tMidPtwMean_AllSk_Area_Cell = tMidPtwMean_SkPxMean_fArea_AllSk_Cell;

% = Vessel Lumen Radius (px), Diameter (px), Area (px2)
oRadius.L.tMidPt_SkPx_RadiusL = tMidPt_SkPx_lRadiusL;
oRadius.L.tMidPt_SkPx_Diameter = tMidPt_SkPx_lDiameter;
oRadius.L.tMidPtwMean_SkPx_Diameter = tMidPtwMean_SkPx_lDiameter;
oRadius.L.tMidPtMin_SkPx_Diameter = tMidPtMin_SkPx_lDiameter;
oRadius.L.tMidPtMax_SkPx_Diameter = tMidPtMax_SkPx_lDiameter;

oRadius.L.tMidPt_Diameter_Cell = tMidPt_lDiameter_Cell;
oRadius.L.tMidPt_AllSk_Diameter_Cell = tMidPt_SkPx_lDiameter_AllSk_Cell;
oRadius.L.tMidPtwMean_SkPxMean_AllSk_Diameter_Cell = tMidPtwMean_SkPxMean_lDiameter_AllSk_Cell;
oRadius.L.tMidPtwMean_SkPxMin_AllSk_Diameter_Cell = tMidPtwMean_SkPxMin_lDiameter_AllSk_Cell;
oRadius.L.tMidPtwMean_SkPxMax_AllSk_Diameter_Cell = tMidPtwMean_SkPxMax_lDiameter_AllSk_Cell;
oRadius.L.tMidPtMax_SkPxMean_AllSk_Diameter_Cell = tMidPtMax_SkPxMean_lDiameter_AllSk_Cell;
oRadius.L.tMidPtMax_SkPxMin_AllSk_Diameter_Cell = tMidPtMax_SkPxMin_lDiameter_AllSk_Cell;
oRadius.L.tMidPtMax_SkPxMax_AllSk_Diameter_Cell = tMidPtMax_SkPxMax_lDiameter_AllSk_Cell;

oRadius.L.tMidPt_Area_Cell = tMidPt_lArea_Cell;
oRadius.L.tMidPt_AllSk_Area_Cell = tMidPt_lArea_AllSk_Cell;
oRadius.L.tMidPtwMean_AllSk_Area_Cell = tMidPtwMean_SkPxMean_lArea_AllSk_Cell;
oRadius.L.tMidPtMax_AllSk_Area_Cell = tMidPtMax_SkPxMean_lArea_AllSk_Cell;


%% Create "skeletons" of even spacing between original skeleton and Img_Thresh Perimeter

[Skel_Crop_xy_Cell,Skel_Crop_bdDist_Cell,Skel_Crop_xy_Flag_Cell] = ...
    fCreateSkelfromThreshPerim(ThreshPerim_Crop_xy_tpCell,ThreshPerim_Crop_xy_Flag_tpCell,tMidPt_SkPx_fRadiusL,NumSkel_VsWidth,0.2*CellDiameter_px);


%% Create Img_Skel_BlkIdx_Crop_Cell w/ new Skeleton Coordinate System

Img_Skel_BlkIdx_Crop_Cell = cell(NumSkel_VsWidth,1,1);
rp_Skel_BlkIdx_Crop_Cell = cell(NumSkel_VsWidth,1,1);

for SkelCt = 1:NumSkel_VsWidth
    Img_Skel_BlkIdx_Crop_Cell{SkelCt,1,1} = zeros(size(ImgStack_Crop(:,:,1)),'uint8');
    for BlkCt = 1:NumBlk
        Img_Skel_BlkIdx_Crop_Cell{SkelCt,1,1}(sub2ind(size(ImgStack_Crop(:,:,1)),Skel_Crop_xy_Cell{SkelCt,1,BlkCt}(:,2),Skel_Crop_xy_Cell{SkelCt,1,BlkCt}(:,1))) = BlkCt; 
    end
    rp_Skel_BlkIdx_Crop_Cell{SkelCt,1,1} = ...
        regionprops(logical(Img_Skel_BlkIdx_Crop_Cell{SkelCt,1,1}),Img_Skel_BlkIdx_Crop_Cell{SkelCt,1,1},'Area','PixelList','PixelIdxList','PixelValues');
end

clearvars Img_Skel_Crop rp_Skel_Crop; 
clearvars SkelCt SegCt MidPtCt BlkCt SkelPxCt;


%% Plot new Skeleton Coordinate System and its relationship w/ Img_Thresh_Crop

figure;
k_Img = 1; 
Img1 = labeloverlay(imadjust(ImgStack_Crop(:,:,k_Img)),bwperim(Img_Thresh_Crop),'Colormap',[1 0 0]);
Img2 = labeloverlay(imadjust(ImgStack_Crop(:,:,k_Img)),max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell,[2,3,1]))),[],3),'Colormap',[0 1 0]);
Img3 = zeros(height(ImgStack_Crop),width(ImgStack_Crop),3);
Img3(:,:,1) = 255.*bwperim(Img_Thresh_Crop);
Img3(:,:,2) = 255.*max(single(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell,[2,3,1]))),[],3)./NumBlk;
Img3(:,:,3) = 255.*bwperim(Img_ThreshVL_Crop);
Img4 = imadjust(ImgStack_Crop(:,:,k_Img));
Img4(logical(rgb2gray(Img3))) = 0;
Img4 = uint8(255.*repmat(Img4,1,1,3))+uint8(1.0*Img3);
Img5 = repmat(255.*uint8(max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell(1,1),[2,3,1]))),[],3)),1,1,3); 
Img6 = repmat(255.*uint8(max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell,[2,3,1]))),[],3)),1,1,3); 
Img7 = repmat(255.*uint8(max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell(1:NumSkel,1),[2,3,1]))),[],3)),1,1,3); 
Img8 = repmat(255.*uint8(bwperim(fConvertSkelxytoMask(logical(max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell,[2,3,1]))),[],3)),find(logical(max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell,[2,3,1]))),[],3))),'n'))),1,1,3); 
Img = horzcat(repmat(255.*imadjust(ImgStack_Crop(:,:,k_Img)),1,1,3),Img1,Img2,Img3,Img4,Img5,Img6,Img7,Img8);
imshow(Img,'InitialMagnification',300);
imwrite(Img,strcat(SaveMODFilePath,'SkeletonCoordinateSystem_RefSlice',num2str(k_Img,'%05d'),'_',timestamp,'.tif'));

% = For Demo Illustrations only
if DEMOOption == 'y'

    Img_Illus_PreROICrop_SCS = zeros(oROI.ImgStack_PreROICrop_Sz(1),oROI.ImgStack_PreROICrop_Sz(2),'logical'); 
    Img_Illus_PreROICrop_SCS(oROI.ROICrop_cuboid(2):oROI.ROICrop_cuboid(2)+height(ImgStack_Crop)-1,...
        oROI.ROICrop_cuboid(1):oROI.ROICrop_cuboid(1)+width(ImgStack_Crop)-1) = max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell,[2,3,1]))),[],3);
    imwrite(255.*uint8(Img_Illus_PreROICrop_SCS),strcat(SaveMODDEMOFilePath,'Img_Illus_PreROI_Crop_SCS_',timestamp,'.tif'));

    Img_Illus_PreROICrop_SCSv = zeros(oROI.ImgStack_PreROICrop_Sz(1),oROI.ImgStack_PreROICrop_Sz(2),'logical');
    Img_Illus_PreROICrop_SCSv(oROI.ROICrop_cuboid(2):oROI.ROICrop_cuboid(2)+height(ImgStack_Crop)-1,...
        oROI.ROICrop_cuboid(1):oROI.ROICrop_cuboid(1)+width(ImgStack_Crop)-1) = max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell(1:NumSkel,1),[2,3,1]))),[],3);
    imwrite(255.*uint8(Img_Illus_PreROICrop_SCSv),strcat(SaveMODDEMOFilePath,'Img_Illus_PreROI_Crop_SCSv_',timestamp,'.tif'));

    Img_Illus_PreROICrop_SCSbyBlk = zeros(oROI.ImgStack_PreROICrop_Sz(1),oROI.ImgStack_PreROICrop_Sz(2),'single');
    Img_Illus_PreROICrop_SCSbyBlk(oROI.ROICrop_cuboid(2):oROI.ROICrop_cuboid(2)+height(ImgStack_Crop)-1,...
        oROI.ROICrop_cuboid(1):oROI.ROICrop_cuboid(1)+width(ImgStack_Crop)-1) = max(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell,[2,3,1])),[],3);
    imwrite(uint8(round(255.*(Img_Illus_PreROICrop_SCSbyBlk./max(Img_Illus_PreROICrop_SCSbyBlk,[],"all")))),strcat(SaveMODDEMOFilePath,'Img_Illus_PreROI_Crop_SCSbyBlk_',timestamp,'.tif'));

    Img_Illus_PreROICrop_SCSPerim = zeros(oROI.ImgStack_PreROICrop_Sz(1),oROI.ImgStack_PreROICrop_Sz(2),'logical');
    Img_Illus_PreROICrop_SCSPerim(oROI.ROICrop_cuboid(2):oROI.ROICrop_cuboid(2)+height(ImgStack_Crop)-1,...
        oROI.ROICrop_cuboid(1):oROI.ROICrop_cuboid(1)+width(ImgStack_Crop)-1) = logical(max(logical(cell2mat(permute(Img_Skel_BlkIdx_Crop_Cell,[2,3,1]))),[],3));
    [Img_Illus_PreROICrop_SCSPerim,~] = fConvertSkelxytoMask(Img_Illus_PreROICrop_SCSPerim,find(Img_Illus_PreROICrop_SCSPerim),'n');
    imwrite(255.*uint8(Img_Illus_PreROICrop_SCSPerim),strcat(SaveMODDEMOFilePath,'Img_Illus_PreROI_Crop_SCSMask_',timestamp,'.tif'));
    imwrite(255.*uint8(bwperim(Img_Illus_PreROICrop_SCSPerim)),strcat(SaveMODDEMOFilePath,'Img_Illus_PreROI_Crop_SCSPerim_',timestamp,'.tif'));

    Img_Illus_PreROICrop_Skel = zeros(oROI.ImgStack_PreROICrop_Sz(1),oROI.ImgStack_PreROICrop_Sz(2),'logical');
    Img_Illus_PreROICrop_Skel(oROI.ROICrop_cuboid(2):oROI.ROICrop_cuboid(2)+height(ImgStack_Crop)-1,...
        oROI.ROICrop_cuboid(1):oROI.ROICrop_cuboid(1)+width(ImgStack_Crop)-1) = logical(Img_Skel_BlkIdx_Crop_Cell{1,1});
    imwrite(255.*uint8(Img_Illus_PreROICrop_Skel),strcat(SaveMODDEMOFilePath,'Img_Illus_PreROI_Crop_Skel_',timestamp,'.tif'));

end

close all;

clearvars k_Img Img Img1 Img2 Img3 Img4 Img5 Img6 Img7 Img8 Img_Illus_PreROICrop_SCS Img_Illus_PreROICrop_SCSv Img_Illus_PreROICrop_SCSbyBlk Img_Illus_PreROICrop_SCSPerim Img_Illus_PreROICrop_Skel SaveMODDEMOFilePath DEMOOption;


%% Create Non-Flow Area and Tissue Masks

[Img_ThreshROI_Crop,ImgStack_ThreshROI_Crop] = ...
    fCreateThreshROIMask(Img_Thresh_Crop_Cell,oROI.Img_Mask,TimeParam);

[Img_ThreshVLROI_Crop,ImgStack_ThreshVLROI_Crop] = ...
    fCreateThreshROIMask(Img_ThreshVL_Crop_Cell,oROI.Img_Mask,TimeParam);

clearvars Img_Thresh_Crop_Cell Img_ThreshVL_Crop_Cell;


%% Create ROI Contrast Map

fprintf('Creating ROI Contrast Map...\n');

% = Size Map
[Img_VesselDiameter_Idx,Img_VesselDiameter_px,Img_VesselDiameter_Idx_Illus] = ...
    fCreateImgStackVesselDiameterMap(ImgStack_Crop,Img_ThreshROI_Crop,Skel_Crop_xy_Cell,tMidPt_SkPx_lDiameter,XY_PxLength);

% = Focus Map
[ImgStack_Crop,Img_Focus_Idx,Img_Focus_Idx_Illus] = ...
    fCreateImgStackFocusQualityMap(ImgStack_Crop,logical(Img_VesselDiameter_px).*Img_ThreshROI_Crop,XY_PxLength);

% = Contrast Map
Img_FocusDiameter_Idx = (Img_Focus_Idx+(1/3)*Img_VesselDiameter_Idx)/(4/3);

% = Save montage of mask
imwrite(uint8(Img_Focus_Idx_Illus),strcat(SaveMODFilePath,'Img_FocusIdx_Illus_',timestamp,'.tif'));
imwrite(double(Img_VesselDiameter_Idx_Illus),strcat(SaveMODFilePath,'Img_VesselDiameterIdx_Illus_',timestamp,'.tif'));
imwrite(double(ind2rgb(round(Img_FocusDiameter_Idx/10.*255), turbo(256))),strcat(SaveMODFilePath,'Img_FocusDiameter_Illus_',timestamp,'.tif')); 

% =  Rank Vessel Block by FocusQuality measured by FocusDiameter_Idx
[Skel_LinIdx_Mtx,~,~] = fLinearizeSkelxyCell(ImgStack_Crop,Skel_Crop_xy_Cell);
FocusDiameterIdx_Mean_Blk = zeros(1,1,NumBlk,'single');
for BlkCt = 1:NumBlk
    FocusDiameterIdx_Mean_Blk(1,1,BlkCt) = ...
        mean(Img_FocusDiameter_Idx(Skel_LinIdx_Mtx(SkelBlockLength*(BlkCt-1)+1:SkelBlockLength*(BlkCt),1:ceil(0.5*NumSkel_VsWidth))),'all','omitmissing');
end
[~,FocusDiameterIdx_Mean_BlkRank] = sort(FocusDiameterIdx_Mean_Blk,3,'descend');

clearvars Skel_LinIdx_Mtx BlkCt FocusDiameterIdx_Mean_Blk;


%% Generate Flag Cells: Vessel_FlowRadius_TimeSeg_Flag_Cell, Vessel_FlowRadius_Skel_Flag_Cell

if Vessel_Diameter_TimeDependence == 'n' % No Flag

    tFr_fRadius_1BlkAllSk_Flag_Cell = repmat({repmat(logical(0),TimeParam.tSeg_FrEnd(end),1)},1,1,NumBlk); 
    tMidPt_fRadius_1BlkAllSk_Flag_Cell = repmat({logical(0)},1,NumMidPt,NumBlk); 
    tMidPt_fRadius_Flag_Cell = repmat({logical(0)},NumSkel,NumMidPt,NumBlk);  

    tMidPt_fRadius_1BlkAllSk_Flag_Ct = repmat({logical(0)},1,1,NumBlk); 
    tMidPt_fRadius_1BlkAllSk_Flag_Pct = repmat({logical(0)},1,1,NumBlk); 

    tFr_fRadius_1BlkAllSk_Flag_Ct = repmat({logical(0)},1,1,NumBlk); 
    tFr_fRadius_1BlkAllSk_Flag_Pct = repmat({logical(0)},1,1,NumBlk); 

    tFrMidPt_fRadius_1BlkAllSk_Flag_Ct = repmat({logical(0)},1,1,NumBlk); 
    tFrMidPt_fRadius_1BlkAllSk_Flag_Pct = repmat({logical(0)},1,1,NumBlk);

    tMidPt_fRadius_AllBlkAllSk_Flag_Ct = 0; 
    tMidPt_fRadius_AllBlkAllSk_Flag_Pct = 0;

    tFr_fRadius_AllBlkAllSk_Flag_Ct = 0;
    tFr_fRadius_AllBlkAllSk_Flag_Pct = 0; 

    tFrMidPt_fRadius_AllBlkAllSk_Flag_Ct = 0; 
    tFrMidPt_fRadius_AllBlkAllSk_Flag_Pct = 0;

else % [Removed]

end

clearvars fRadiusL_Frac_FlagID;


%% Generate Flag Cells: Vessel_Radius_TimeSeg_Flag_Cell, Vessel_Radius_Skel_Flag_Cell

if Vessel_Diameter_TimeDependence == 'n' % Dummy

    tFr_lRadius_1BlkAllSk_Flag_Cell = repmat({repmat(logical(0),TimeParam.tSeg_FrEnd(end),1)},1,1,NumBlk); 
    tMidPt_lRadius_1BlkAllSk_Flag_Cell = repmat({logical(0)},1,NumMidPt,NumBlk); 
    tMidPt_lRadius_Flag_Cell = repmat({logical(0)},NumSkel,NumMidPt,NumBlk);
    
    tMidPt_lRadius_1BlkAllSk_Flag_Ct = repmat({logical(0)},1,1,NumBlk); 
    tMidPt_lRadius_1BlkAllSk_Flag_Pct = repmat({logical(0)},1,1,NumBlk);

    tFr_lRadius_1BlkAllSk_Flag_Ct = repmat({logical(0)},1,1,NumBlk);
    tFr_lRadius_1BlkAllSk_Flag_Pct = repmat({logical(0)},1,1,NumBlk); 

    tFrMidPt_lRadius_1BlkAllSk_Flag_Ct = repmat({logical(0)},1,1,NumBlk); 
    tFrMidPt_lRadius_1BlkAllSk_Flag_Pct = repmat({logical(0)},1,1,NumBlk); 

    tMidPt_lRadius_AllBlkAllSk_Flag_Ct = 0; 
    tMidPt_lRadius_AllBlkAllSk_Flag_Pct = 0;

    tFr_lRadius_AllBlkAllSk_Flag_Ct = 0; 
    tFr_lRadius_AllBlkAllSk_Flag_Pct = 0;

    tFrMidPt_lRadius_AllBlkAllSk_Flag_Ct = 0;
    tFrMidPt_lRadius_AllBlkAllSk_Flag_Pct = 0;

else % [Removed]

end

clearvars lRadiusL_Frac_FlagID;


%% Store Flag Cells: Vessel_Radius_TimeSeg_Flag_Cell, Vessel_FlowRadius_TimeSeg_Flag_Cell

oFlag_Radius = struct;

% = Flow (f) Profile Radius
oFlag_Radius.F.tFr_1BlkAllSk_Flag_Cell = tFr_fRadius_1BlkAllSk_Flag_Cell; 
oFlag_Radius.F.tMidPt_1BlkAllSk_Flag_Cell = tMidPt_fRadius_1BlkAllSk_Flag_Cell; 
oFlag_Radius.F.tMidPt_Flag_Cell = tMidPt_fRadius_Flag_Cell; 

oFlag_Radius.F.tMidPt_1BlkAllSk_Flag_Ct = tMidPt_fRadius_1BlkAllSk_Flag_Ct;
oFlag_Radius.F.tMidPt_1BlkAllSk_Flag_Pct = tMidPt_fRadius_1BlkAllSk_Flag_Pct;
oFlag_Radius.F.tFr_1BlkAllSk_Flag_Ct = tFr_fRadius_1BlkAllSk_Flag_Ct;
oFlag_Radius.F.tFr_1BlkAllSk_Flag_Pct = tFr_fRadius_1BlkAllSk_Flag_Pct;
oFlag_Radius.F.tFrMidPt_1BlkAllSk_Flag_Ct = tFrMidPt_fRadius_1BlkAllSk_Flag_Ct;
oFlag_Radius.F.tFrMidPt_1BlkAllSk_Flag_Pct = tFrMidPt_fRadius_1BlkAllSk_Flag_Pct;

oFlag_Radius.F.tMidPt_AllBlkAllSk_Flag_Ct = tMidPt_fRadius_AllBlkAllSk_Flag_Ct;
oFlag_Radius.F.tMidPt_AllBlkAllSk_Flag_Pct = tMidPt_fRadius_AllBlkAllSk_Flag_Pct;
oFlag_Radius.F.tFr_AllBlkAllSk_Flag_Ct = tFr_fRadius_AllBlkAllSk_Flag_Ct;
oFlag_Radius.F.tFr_AllBlkAllSk_Flag_Pct = tFr_fRadius_AllBlkAllSk_Flag_Pct;
oFlag_Radius.F.tFrMidPt_AllBlkAllSk_Flag_Ct = tFrMidPt_fRadius_AllBlkAllSk_Flag_Ct;
oFlag_Radius.F.tFrMidPt_AllBlkAllSk_Flag_Pct = tFrMidPt_fRadius_AllBlkAllSk_Flag_Pct;

% = Lumen (L) Radius
oFlag_Radius.L.tFr_1BlkAllSk_Flag_Cell = tFr_lRadius_1BlkAllSk_Flag_Cell; 
oFlag_Radius.L.tMidPt_1BlkAllSk_Flag_Cell = tMidPt_lRadius_1BlkAllSk_Flag_Cell;
oFlag_Radius.L.tMidPt_Flag_Cell = tMidPt_lRadius_Flag_Cell;

oFlag_Radius.L.tMidPt_1BlkAllSk_Flag_Ct = tMidPt_lRadius_1BlkAllSk_Flag_Ct;
oFlag_Radius.L.tMidPt_1BlkAllSk_Flag_Pct = tMidPt_lRadius_1BlkAllSk_Flag_Pct;
oFlag_Radius.L.tFr_1BlkAllSk_Flag_Ct = tFr_lRadius_1BlkAllSk_Flag_Ct;
oFlag_Radius.L.tFr_1BlkAllSk_Flag_Pct = tFr_lRadius_1BlkAllSk_Flag_Pct;
oFlag_Radius.L.tFrMidPt_1BlkAllSk_Flag_Ct = tFrMidPt_lRadius_1BlkAllSk_Flag_Ct;
oFlag_Radius.L.tFrMidPt_1BlkAllSk_Flag_Pct = tFrMidPt_lRadius_1BlkAllSk_Flag_Pct;

oFlag_Radius.L.tMidPt_AllBlkAllSk_Flag_Ct = tMidPt_lRadius_AllBlkAllSk_Flag_Ct;
oFlag_Radius.L.tMidPt_AllBlkAllSk_Flag_Pct = tMidPt_lRadius_AllBlkAllSk_Flag_Pct;
oFlag_Radius.L.tFr_AllBlkAllSk_Flag_Ct = tFr_lRadius_AllBlkAllSk_Flag_Ct;
oFlag_Radius.L.tFr_AllBlkAllSk_Flag_Pct = tFr_lRadius_AllBlkAllSk_Flag_Pct;
oFlag_Radius.L.tFrMidPt_AllBlkAllSk_Flag_Ct = tFrMidPt_lRadius_AllBlkAllSk_Flag_Ct;
oFlag_Radius.L.tFrMidPt_AllBlkAllSk_Flag_Pct = tFrMidPt_lRadius_AllBlkAllSk_Flag_Pct;

clearvars fRadius_Flag_Cell lRadius_Flag_Cell;
clearvars or_wrapper;


%% Collect to-be-saved Skeleton Coordinate System (SCS) Parameters in structure

oSCS = struct;

oSCS.SkelBlockLength = SkelBlockLength;

oSCS.NumSkel_VsWidth = NumSkel_VsWidth;
oSCS.NumSkel = NumSkel;
oSCS.NumSeg = NumSeg;
oSCS.NumMidPt= NumMidPt;
oSCS.NumBlk= NumBlk;

oSCS.Skel_xy_Cell = Skel_Crop_xy_Cell; 
oSCS.Skel_xy_Flag_Cell = Skel_Crop_xy_Flag_Cell;  

oSCS.ImgStack_fThreshROI = ImgStack_ThreshROI_Crop; % Time-dependent flow profile ThreshROI
oSCS.Img_fThreshROI = Img_ThreshROI_Crop; % Time-independent flow profile ThreshROI

oSCS.ImgStack_lThreshROI = ImgStack_ThreshVLROI_Crop; % Time-dependent vessel lumen ThreshROI
oSCS.Img_lThreshROI = Img_ThreshVLROI_Crop; % Time-independent vessel lumen ThreshROI

oSCS.Skel_bdDist_Cell = Skel_Crop_bdDist_Cell;

oSCS.Img_Skel_BlkIdx_Cell =  Img_Skel_BlkIdx_Crop_Cell; 
oSCS.rp_Skel_BlkIdx_Cell = rp_Skel_BlkIdx_Crop_Cell;

clearvars Skel_Crop_xy_Cell Skel_Crop_xy_Flag_Cell;
clearvars ImgStack_ThreshROI_Crop Img_ThreshROI_Crop ImgStack_ThreshVLROI_Crop Img_ThreshVLROI_Crop;
clearvars Skel_Crop_bdDist_Cell;
clearvars Img_Skel_BlkIdx_Crop_Cell rp_Skel_BlkIdx_Crop_Cell;


%% Collect to-be-saved Image Contrast Map Parameters in structure

oFocusQ = struct;

oFocusQ.Img_FocusDiameter_Idx = Img_FocusDiameter_Idx; % Final contrast map based on focus and size maps
oFocusQ.Img_Focus_Idx = Img_Focus_Idx; % Focus component of final contrast map
oFocusQ.Img_Diameter_Idx = Img_VesselDiameter_Idx; % Vessel diameter component of final contrast map

oFocusQ.BlkRank = FocusDiameterIdx_Mean_BlkRank;

clearvars Img_FocusDiameter_Idx
clearvars Img_Focus_Idx;
clearvars Img_VesselDiameter_Idx;
clearvars FocusDiameterIdx_Mean_BlkRank;


%% Clear variables likely no longer relevant after this module

clearvars Img_Skel_Crop_sCell Img_Thresh_Crop_sCell; % Obsolete
clearvars rp_Skel_Crop_sCell rp_Thresh_Crop_sCell;

clearvars ThreshPerim_Crop_xy_tpCell ThreshPerim_Crop_xy_Flag_tpCell;
clearvars rp_ThreshPerim_Crop_tpCell NumThreshPerim;

clearvars ThreshVLPerim_Crop_xy_tpCell ThreshVLPerim_Crop_xy_Flag_tpCell;
clearvars rp_ThreshVLPerim_Crop_tpCell NumThreshPerim;

clearvars NumThreshVLPerim;

clearvars Img_Skel_Crop rp_Skel_Crop; 
clearvars Img_Thresh_Crop rp_Thresh_Crop; 
clearvars Img_ROI_Crop rp_ROI_Crop; 

clearvars Img_SkelVL_Crop rp_SkelVL_Crop; 
clearvars Img_ThreshVL_Crop Img_ROI_Crop rp_ThreshVL_Crop rp_ROI_Crop; 

clearvars GraphCut_ROI_Stack_Crop;

clearvars PolygonROI_startXY;

clearvars Img_*_Illus Img_VesselDiameter_px;

clearvars *Radius* -except oRadius oFlag_Radius;
clearvars *Diameter* -except Vessel_Diameter_TimeDependence NumCellDiameter_in_SkelBlockLength;
clearvars *Area*;


%% Save Workspace 

% = OPTION 2: Save only new variables introduced
if ~exist('varList_MOD01') && ~exist('varList_MOD02') 
    varList_MOD01 = who;
    varList_MOD01 = varList_MOD01(~ismember(varList_MOD01,varList_All));
    varList_All = vertcat(varList_All,varList_MOD01);
    save(strcat(SaveMODFilePath,'Workspace_STEP02Var_',SampleIDString,'_',timestamp,'.mat'),varList_MOD01{:},...
        'TimeParam'); 
end


%% Save script in directory

ScriptName=mfilename;
PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveMODFilePath);
publish(strcat(ScriptName,'.m'),PublishOptions);
