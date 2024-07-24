function [Img_BW_Skel_Output,rp_Img_BW_Skel_Output] = fCreateImgSkel(Img_BW_Input,MinSkelBranchLength_Input)

Img_BW_Skel_Output = ... 
    bwskel(Img_BW_Input,'MinBranchLength',MinSkelBranchLength_Input);

% = Remove short skeleton segments
Img_BW_Skel_BranchPt = bwmorph(Img_BW_Skel_Output, 'branchpoints');

% figure;
% imshowpair(Img_Skel,Img_Skel_BranchPt,'falseColor');

Img_BW_Skel_Output(Img_BW_Skel_BranchPt) = 0;
rp_BW_Skel_Output = regionprops(Img_BW_Skel_Output,'Area','PixelList','PixelIdxList','Orientation');

rp_BW_Skel_Output = rp_BW_Skel_Output(cat(1,rp_BW_Skel_Output.Area)>MinSkelBranchLength_Input);

% = Remove short branches. Update Skeleton Image
Img_BW_Skel_Output = zeros(size(Img_BW_Skel_Output),'logical');
Img_BW_Skel_Output(cat(1,rp_BW_Skel_Output.PixelIdxList)) = 1;
Img_BW_Skel_Output = bwmorph(Img_BW_Skel_Output,'bridge'); 
Img_BW_Skel_Output = bwmorph(Img_BW_Skel_Output,'thin'); 

% = Smooth out minor kinks in skeleton
Img_BW_Skel_Output = bwskel(bwmorph(imdilate(Img_BW_Skel_Output,strel('disk',3)),'majority')); 
Img_BW_Skel_Output = bwmorph(Img_BW_Skel_Output,'thin');

rp_BW_Skel_Output = regionprops(Img_BW_Skel_Output,'Area','PixelList','PixelIdxList','Orientation');