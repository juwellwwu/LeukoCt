function [Skel_xy_Ord_Cell_Output,Skel_xy_Flag_Cell_Output,ImgStack_Skel_Output,rp_Skel_Cell_Output] = ...
    fMatchSkelLength(Skel_xy_Ord_Cell_Input,Intersect_DilateSkelPxIdx_SkCell_Input,ImgStack_Skel_Input,rp_Skel_Cell_Input)

% = Fill in NaN @ High Curvature Points

Intersect_DilateSkelPxIdx = zeros(size(Intersect_DilateSkelPxIdx_SkCell_Input),'single'); 
Intersect_DilateSkelPxIdx_NaN = zeros(size(Intersect_DilateSkelPxIdx_SkCell_Input),'logical');

for hCt = 1:height(Intersect_DilateSkelPxIdx_SkCell_Input) 

    % Locate NaN pixels in ThreshPerim designated by hCt
    rp_Temp = regionprops(isnan(cell2mat(Intersect_DilateSkelPxIdx_SkCell_Input(hCt,:))),'PixelIdxList');
    
    % Blobs to exclude from missing
    FillMissing_SkipID = zeros(numel(rp_Temp),1,'logical');
    for bCt = 1:numel(rp_Temp) % # blobs = # NaN regions
        FillMissing_SkipID(bCt,1) = ismember(1,rp_Temp(bCt).PixelIdxList) || ismember(width(Intersect_DilateSkelPxIdx_SkCell_Input),rp_Temp(bCt).PixelIdxList);
    end
    
    % Fill all NaN values in ThreshPerim with linear interpolation of nearest neighbor
    Intersect_DilateSkelPxIdx(hCt,:) = fillmissing(cell2mat(Intersect_DilateSkelPxIdx_SkCell_Input(hCt,:)),'linear',2);
    
    % Fill back all pixels flagged by FillMissing_SkipID (at end points) w/ NaN
    Intersect_DilateSkelPxIdx(hCt,cat(1,rp_Temp(FillMissing_SkipID).PixelIdxList)) = NaN;
    
    % Record where the pixels were NaN in Intersect_DilateSkelPxIdx_SkCell_Input
    Intersect_DilateSkelPxIdx_NaN(hCt,:) = isnan(cell2mat(Intersect_DilateSkelPxIdx_SkCell_Input(hCt,:)));

end

% Remove ThreshPerimPxCt (Col) at two endpoints with NaN @ any ThreshPerim
% (Row)
Intersect_DilateSkelPxIdx_NaN(:,isnan(max(Intersect_DilateSkelPxIdx,[],1,"includenan"))) = [];
Intersect_DilateSkelPxIdx(:,isnan(max(Intersect_DilateSkelPxIdx,[],1,"includenan"))) = [];

% = Update ThreshPerims
SkelPxIdx_Cell = cell(height(Skel_xy_Ord_Cell_Input),1); 
Img_Temp_Cell = cell(height(Skel_xy_Ord_Cell_Input),1); 

Skel_xy_Flag_Cell_Output = cell(height(Skel_xy_Ord_Cell_Input),1); 

ImgStack_Skel_Output = ImgStack_Skel_Input;
rp_Skel_Cell_Output = rp_Skel_Cell_Input;
Skel_xy_Ord_Cell_Output = Skel_xy_Ord_Cell_Input;

for SkelCt = 1:height(Skel_xy_Ord_Cell_Input)

    SkelPxIdx_Cell{SkelCt,1} = round(Intersect_DilateSkelPxIdx(SkelCt,:));
    if ~all((SkelPxIdx_Cell{SkelCt,1})>0) % Debug
        pause()
    end
    Skel_xy_Ord_Cell_Output{SkelCt,1} = Skel_xy_Ord_Cell_Output{SkelCt,1}(SkelPxIdx_Cell{SkelCt,1},:); 

    Skel_xy_Flag_Cell_Output{SkelCt,1} = Intersect_DilateSkelPxIdx_NaN(SkelCt,:)';

    Img_Temp_Cell{SkelCt,1} = zeros(size(ImgStack_Skel_Output(:,:,SkelCt)),'logical');
    Img_Temp_Cell{SkelCt,1}(sub2ind(size(ImgStack_Skel_Output(:,:,SkelCt)),Skel_xy_Ord_Cell_Input{SkelCt,1}(:,2),Skel_xy_Ord_Cell_Input{SkelCt,1}(:,1))) = 1;

    ImgStack_Skel_Output(:,:,SkelCt) = Img_Temp_Cell{SkelCt,1};
    rp_Skel_Cell_Output{SkelCt,1} = regionprops(ImgStack_Skel_Output(:,:,SkelCt),'Area','PixelList','PixelIdxList');

end