function [Img_BW_SmoothBoundary_Output,Img_BW_Smooth_Output] = fSmoothImgBWBoundary_SG(Img_BW_Input,SGFilter_Order,SGFilter_FrLength,XY_PxLength_Input)

CellDiameter_px_Input = 5/XY_PxLength_Input;

Img_Boundary = bwboundaries(Img_BW_Input);

Img_Boundary_Length_px = cell2mat(cellfun(@height,Img_Boundary,'UniformOutput',false));
Img_Boundary = Img_Boundary(Img_Boundary_Length_px>100);

% Set Savitzky-Golay sliding polynomial filter parameters:
sgolayfilt_X_wrapper = @(x) sgolayfilt(x(:,2),SGFilter_Order,SGFilter_FrLength);
sgolayfilt_Y_wrapper = @(x) sgolayfilt(x(:,1),SGFilter_Order,SGFilter_FrLength);

Img_Boundary_smoothX = cellfun(sgolayfilt_X_wrapper,Img_Boundary,'UniformOutput', false);
Img_Boundary_smoothY = cellfun(sgolayfilt_Y_wrapper,Img_Boundary,'UniformOutput', false);

% Combine boundaries 
Img_Boundary_smoothX = cell2mat(Img_Boundary_smoothX);
Img_Boundary_smoothY = cell2mat(Img_Boundary_smoothY);

% Remove out of bounds subscript
Img_Boundary_smoothY((Img_Boundary_smoothX<0.5)|(Img_Boundary_smoothX>(size(Img_BW_Input,2)+0.5-eps))) = [];
Img_Boundary_smoothX((Img_Boundary_smoothX<0.5)|(Img_Boundary_smoothX>(size(Img_BW_Input,2)+0.5-eps))) = [];

Img_Boundary_smoothX((Img_Boundary_smoothY<0.5)|(Img_Boundary_smoothY>(size(Img_BW_Input,1)+0.5-eps))) = [];
Img_Boundary_smoothY((Img_Boundary_smoothY<0.5)|(Img_Boundary_smoothY>(size(Img_BW_Input,1)+0.5-eps))) = [];

% Convert smoothed X and Y to linear index
Img_Boundary_smoothLinIdx = ...
    sub2ind(size(Img_BW_Input),round(Img_Boundary_smoothY),round(Img_Boundary_smoothX));

% Create Image
Img_BW_SmoothBoundary_Output = zeros(size(Img_BW_Input),'logical');
Img_BW_SmoothBoundary_Output(Img_Boundary_smoothLinIdx) = 1;
Img_BW_SmoothBoundary_Output = imdilate(Img_BW_SmoothBoundary_Output,strel('disk',1)); 

% Create smoothed mask
Img_BW_imfillSkel = bwmorph(bwskel(Img_BW_Input),'spur',3);
Img_BW_imfillSkel = Img_BW_imfillSkel.*~Img_BW_SmoothBoundary_Output;
rp_BW_imfillSkel = regionprops(Img_BW_imfillSkel,'Area','PixelIdxList');
FlagID = cat(1,rp_BW_imfillSkel.Area) < 4*CellDiameter_px_Input;
Img_BW_imfillSkel(cat(1,rp_BW_imfillSkel.Area)) = 0;
Img_BW_Smooth_Output = ...
    logical(imfill(Img_BW_SmoothBoundary_Output,find(Img_BW_imfillSkel))); 
