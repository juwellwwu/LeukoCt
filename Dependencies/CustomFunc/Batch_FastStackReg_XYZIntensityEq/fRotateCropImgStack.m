function [ImgStack_RotCrop] = fRotateCropImgStack(ImgStack_In, deg_Rot_In)

%% Create mask that specifies where 0 values from rotation is 

Img_Rot_Mask_PreClean = ...
    logical(imrotate(ones(size(ImgStack_In(:,:,1)),'logical'),deg_Rot_In,"nearest","loose"));

Img_Rot_Mask = Img_Rot_Mask_PreClean;

Img_Rot_Mask(1:ceil(sind(deg_Rot_In)*width(ImgStack_In)),:) = 0;
Img_Rot_Mask(end-ceil(sind(deg_Rot_In)*width(ImgStack_In))+1:end,:) = 0;

Img_Rot_Mask(:,1:ceil(sind(deg_Rot_In)*height(ImgStack_In))) = 0;
Img_Rot_Mask(:,end-ceil(sind(deg_Rot_In)*height(ImgStack_In))+1:end) = 0;

rp_Rot_Mask = regionprops(Img_Rot_Mask,'BoundingBox');

% Check accuracy of 0 values removal
% figure; 
% imshowpair(Img_Rot_Mask_PreClean,Img_Rot_Mask);


%% Rotate and Crop

% figure;
ImgStack_RotCrop = ...
    zeros(rp_Rot_Mask.BoundingBox(4)+1,rp_Rot_Mask.BoundingBox(3)+1,'single');
for k =1:size(ImgStack_In,3)
    Img_Rot = imrotate(ImgStack_In(:,:,k),deg_Rot_In,"nearest","loose");
    ImgStack_RotCrop(:,:,k) = imcrop(Img_Rot,rp_Rot_Mask.BoundingBox);
    % imshow(ImgStack_RotCrop(:,:,k));
end


