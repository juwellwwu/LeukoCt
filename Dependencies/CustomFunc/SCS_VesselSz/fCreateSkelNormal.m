function [SkelNorm_xlim,SkelNorm_ylim,SkelNorm_FullLn_Endpt_RowIdx,SkelCurvature_px] = ...
    fCreateSkelNormal(Skel_xy_Ord_Input,MovWindowSize,NormHalfLength)

% = Determine moving window indices for skeleton normal and curvature
% calculations
% LineNormals2D() and LineCurvature2D() from FileExchange: 
% https://www.mathworks.com/matlabcentral/fileexchange/32696-2d-line-curvature-and-normals

% Calculate Ln start and end points
SkelNorm_Ln_Start = (-MovWindowSize:1:height(Skel_xy_Ord_Input)); 
SkelNorm_Ln_End = SkelNorm_Ln_Start+MovWindowSize-1;
SkelNorm_Ln_Midpt = round((SkelNorm_Ln_Start+SkelNorm_Ln_End)*0.5); 

SkelNorm_Ln = horzcat(SkelNorm_Ln_Start',SkelNorm_Ln_End',SkelNorm_Ln_Midpt');

SkelNorm_Ln(SkelNorm_Ln(:,1)<1,1)=1;
SkelNorm_Ln(SkelNorm_Ln(:,2)>height(Skel_xy_Ord_Input),2)=height(Skel_xy_Ord_Input);

SkelNorm_Ln(SkelNorm_Ln(:,3)>height(Skel_xy_Ord_Input),:) = [];
SkelNorm_Ln(SkelNorm_Ln(:,3)<1,:) = [];

% = Calculate normals
Skel_Norm = LineNormals2D(Skel_xy_Ord_Input, SkelNorm_Ln(:,1:2)); % [x,y] of skeleton

% = Calculate, for the normals of all skeleton pixels, xlimits and ylimits of the normal line. 
SkelNorm_xlim = ... 
    [Skel_xy_Ord_Input(:,1)-(NormHalfLength*Skel_Norm(:,1)) Skel_xy_Ord_Input(:,1)+(NormHalfLength*Skel_Norm(:,1))];
SkelNorm_ylim = ... 
    [Skel_xy_Ord_Input(:,2)-(NormHalfLength*Skel_Norm(:,2)) Skel_xy_Ord_Input(:,2)+(NormHalfLength*Skel_Norm(:,2))];

% = Find row nearest to skeleton endpoint with full moving window for normals measurement
SkelNorm_FullLn_Endpt_RowIdx = [find(diff(SkelNorm_Ln(:,1)),1,'first'), find(diff(SkelNorm_Ln(:,2)),1,'last')];

% == Update other skeleton line normal parameters
% Cut off endpoints without complete moving window for line normal
SkelNorm_xlim([1:(SkelNorm_FullLn_Endpt_RowIdx(1)-1),(SkelNorm_FullLn_Endpt_RowIdx(2)+1):end],:) = [];
SkelNorm_ylim([1:(SkelNorm_FullLn_Endpt_RowIdx(1)-1),(SkelNorm_FullLn_Endpt_RowIdx(2)+1):end],:) = [];

% = Calculate skeleton curvature
SkelCurvature_px = LineCurvature2D(Skel_xy_Ord_Input, SkelNorm_Ln(:,1:2)); % [x,y] of skeleton
SkelCurvature_px([1:(SkelNorm_FullLn_Endpt_RowIdx(1)-1),(SkelNorm_FullLn_Endpt_RowIdx(2)+1):end],:) = [];

