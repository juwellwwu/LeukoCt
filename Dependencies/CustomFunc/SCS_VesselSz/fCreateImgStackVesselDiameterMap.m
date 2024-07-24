function [Img_VesselDiameter_Idx_Out,Img_VesselDiameter_Out,Img_VesselDiameter_Idx_Illus_Out] = ...
    fCreateImgStackVesselDiameterMap(ImgStack_In,ImgStack_TissueMask_In,Skel_xy_Cell_In,tMidPt_Vessel_Diameter_px_In,XY_PxLength_In)

% Re-organize Skel_xy_Cell_In
[Skel_LinIdx_Mtx_In,Skel_x_Mtx_In,Skel_y_Mtx_In] = fLinearizeSkelxyCell(ImgStack_In,Skel_xy_Cell_In);

% = Create Map
Img_VesselDiameter_Out = nan(height(ImgStack_In),width(ImgStack_In));

% Replace skeleton pixels with diameter value
for SkelPxCt = 1:height(Skel_LinIdx_Mtx_In)
    Img_VesselDiameter_Out(Skel_LinIdx_Mtx_In(SkelPxCt,:)) = max(tMidPt_Vessel_Diameter_px_In(SkelPxCt,:),[],2); 
end

% Fill NaN values by interpolation
Img_VesselDiameter_Out = fillmissing2(Img_VesselDiameter_Out,'linear'); 

% Fill tissue space with 0
Img_VesselDiameter_Out(~ImgStack_TissueMask_In) = 0;

% Set up Vessel Diameter Index image
% Index = 10 for up to Diameter = 15, then -1 every 1.5 um
% Index = 10 from Vessel Diameter <= 15 um
% Index = 9 from Vessel Diameter 15-16.5
% Index = 8 from Vessel Diameter 16.5-18 ...
% Index = 7 from Vessel Diameter 18.0-19.5 ...
% Index = 6 from Vessel Diameter 19.5-21.0 ...
% Index = 5 from Vessel Diameter 21.0-22.5 ...
% Index = 4 from Vessel Diameter 22.5-24.0 ...
% Index = 3 from Vessel Diameter 24.0-25.5 ...
% Index = 2 from Vessel Diameter 25.5-27.0 ...
% Index = 1 from Vessel Diameter 27-28.5 
% Index = 0 from Vessel Diameter >28.5

Img_VesselDiameter_Idx_Out = zeros(size(Img_VesselDiameter_Out),'single');

Diameter_Cutoff = (15:1.5:28.5)/XY_PxLength_In;
Diameter_Idx = (9:-1:1);

Img_VesselDiameter_Idx_Out(Img_VesselDiameter_Out>eps & Img_VesselDiameter_Out<=Diameter_Cutoff(1)) = 10; 
for idxCt = 1:numel(Diameter_Idx)
    Img_VesselDiameter_Idx_Out(Img_VesselDiameter_Out>Diameter_Cutoff(idxCt) & Img_VesselDiameter_Out<=Diameter_Cutoff(idxCt+1)) = Diameter_Idx(idxCt);
end

% Set NaN values = 0;
Img_VesselDiameter_Out(isnan(Img_VesselDiameter_Out)) = 0;

% Check
% figure;
% imshow(Img_VesselDiameter_Idx_Out./max(Img_VesselDiameter_Idx_Out,[],"all"));

Img_VesselDiameter_Idx_Illus_Out = horzcat(Img_VesselDiameter_Out./max(Img_VesselDiameter_Out,[],'all'),Img_VesselDiameter_Idx_Out/10);

