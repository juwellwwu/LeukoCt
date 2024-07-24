function [ImgStack_RegCrop_Output,RegShiftbyZ_Output] = fCropImgStackReg(ImgStack_Reg_Input,RegShiftbyZ_Input,RegShiftbyZScale,XY_PxLength_Input,FlowMaskOption)

fprintf('Crop registered image stack to remove affected areas...\n');

% = Define volume to crop based on RegShiftbyZ 
RegShiftbyZ_Output = RegShiftbyZ_Input; 

NumPx_Top = ceil(abs(max(RegShiftbyZ_Output(RegShiftbyZ_Output(:,3)>0,3)))); 
NumPx_Bttm = ceil(abs(min(RegShiftbyZ_Output(RegShiftbyZ_Output(:,3)<0,3)))); 
NumPx_Left = ceil(abs(max(RegShiftbyZ_Output(RegShiftbyZ_Output(:,4)>0,4)))); 
NumPx_Right = ceil(abs(min(RegShiftbyZ_Output(RegShiftbyZ_Output(:,4)<0,4)))); 

if isempty(NumPx_Top) 
    NumPx_Top = 1;
end

if isempty(NumPx_Bttm)
    NumPx_Bttm = 1;
end

if isempty(NumPx_Left) 
    NumPx_Left = 1;
end

if isempty(NumPx_Right)
    NumPx_Right = 1;
end

cuboid_RegShiftbyZ = ... 
    [NumPx_Left NumPx_Top 1 ... 
    size(ImgStack_Reg_Input,2)-NumPx_Left-NumPx_Right size(ImgStack_Reg_Input,1)-NumPx_Top-NumPx_Bttm size(ImgStack_Reg_Input,3)];

Img_RegShiftbyZ = zeros(height(ImgStack_Reg_Input),width(ImgStack_Reg_Input),'logical');
if ~isequal(RegShiftbyZScale,1)
    Img_RegShiftbyZ = logical(imresize(Img_RegShiftbyZ,RegShiftbyZScale,'Method','bilinear')); 
end
Img_RegShiftbyZ = ... 
        insertShape(uint8(Img_RegShiftbyZ),"filled-rectangle",[cuboid_RegShiftbyZ(1) cuboid_RegShiftbyZ(2) cuboid_RegShiftbyZ(4) cuboid_RegShiftbyZ(5)],'Color','white');
Img_RegShiftbyZ = imbinarize(mean(Img_RegShiftbyZ,3)); 

Img_RegShiftbyZ = imresize(Img_RegShiftbyZ,[height(ImgStack_Reg_Input) width(ImgStack_Reg_Input)],'Method',"bilinear");  

% figure;
% imshow(Img_RegShiftbyZ);

rp_RegShiftbyZ = regionprops(Img_RegShiftbyZ,'BoundingBox','PixelIdxList'); 

if isempty(rp_RegShiftbyZ) 
    ImgStack_RegCrop_Output  = ImgStack_Reg_Input;
    RegShiftbyZ_Output = RegShiftbyZ_Input;
    return
end

cuboid_RegShiftbyZ = ... 
    [ceil(rp_RegShiftbyZ(1).BoundingBox(1)) ceil(rp_RegShiftbyZ(1).BoundingBox(2)) 1 ... 
    floor(rp_RegShiftbyZ(1).BoundingBox(3)) floor(rp_RegShiftbyZ(1).BoundingBox(4)) size(ImgStack_Reg_Input,3)];

rp_RegShiftbyZ = regionprops(imcomplement(Img_RegShiftbyZ),'BoundingBox','PixelIdxList'); 

% = Verify that Crop Box does not touch Flow Areas
if FlowMaskOption == 'y'

    [~,~,Img_Mask_Flow] = fCreateRegImgZthresh(ImgStack_Reg_Input,XY_PxLength_Input,"",'n','n','n');
    
    rp_Mask = regionprops(logical(Img_Mask_Flow),'Area','boundingbox','Centroid','PixelList','PixelIdxList');
    
    FlagID = zeros(height(rp_Mask),1,'logical');
    for BlobCt = 1:height(rp_Mask)
        FlagID(BlobCt,1) = numel(find(ismember(rp_Mask(BlobCt).PixelIdxList,rp_RegShiftbyZ(1).PixelIdxList))) > 0.5*rp_Mask(BlobCt).Area; % ismember(A,B) = A found in B
    end

    % Check
    % % % Img_Mask_Flow(cat(1,rp_Mask(FlagID).PixelIdxList)) = 0;
    % % % figure;
    % % % imshow(Img_Mask_Flow);

    rp_Mask(FlagID) = [];

    rp_Mask_Bbox = cat(1,rp_Mask.BoundingBox);
    rp_Mask_Bbox = ... 
        horzcat(rp_Mask_Bbox,rp_Mask_Bbox(:,1)+rp_Mask_Bbox(:,3),rp_Mask_Bbox(:,2)+rp_Mask_Bbox(:,4));

    rp_Mask_Bbox_Allxmin_Allymin = min(rp_Mask_Bbox(:,1:2),[],1); 
    rp_Mask_Bbox_Allxmin_Allymin = rp_Mask_Bbox_Allxmin_Allymin-30; 
    rp_Mask_Bbox_Allxmin_Allymin(rp_Mask_Bbox_Allxmin_Allymin<1) = 1;

    rp_Mask_Bbox_Allxmax_Allymax = max(rp_Mask_Bbox(:,5:6),[],1); 
    rp_Mask_Bbox_Allxmax_Allymax = rp_Mask_Bbox_Allxmax_Allymax+30; 
    rp_Mask_Bbox_Allxmax_Allymax(1,1) = ... 
        min(rp_Mask_Bbox_Allxmax_Allymax(1,1),width(ImgStack_Reg_Input));
    rp_Mask_Bbox_Allxmax_Allymax(1,2) = ... 
        min(rp_Mask_Bbox_Allxmax_Allymax(1,2),height(ImgStack_Reg_Input));

end

% = Combine two cuboid information

if FlowMaskOption == 'y'

    cuboid_xmin = min(ceil(rp_RegShiftbyZ(1).BoundingBox(1)),rp_Mask_Bbox_Allxmin_Allymin(1,1));
    cuboid_ymin = min(ceil(rp_RegShiftbyZ(1).BoundingBox(2)),rp_Mask_Bbox_Allxmin_Allymin(1,2));
    cuboid_zmin = 1;

    cuboid_xmax = max(ceil(rp_RegShiftbyZ(1).BoundingBox(1))+floor(rp_RegShiftbyZ(1).BoundingBox(3))-1,rp_Mask_Bbox_Allxmax_Allymax(1,1));
    cuboid_ymax = max(ceil(rp_RegShiftbyZ(1).BoundingBox(2))+floor(rp_RegShiftbyZ(1).BoundingBox(4))-1,rp_Mask_Bbox_Allxmax_Allymax(1,2));
    cuboid_zmax = size(ImgStack_Reg_Input,3);

    cuboid = [cuboid_xmin cuboid_xmin cuboid_zmin ...
    cuboid_xmax-cuboid_xmin cuboid_ymax-cuboid_ymin cuboid_zmax-cuboid_zmin];

else

    cuboid = cuboid_RegShiftbyZ;

end

% = Finally: ensure size(ImStack_Reg_Input) not violated
% MATLAB treats each pixel pt as spanning (Px-0.5,Px+0.5) in x,y
if (cuboid(1)+cuboid(4)) > (width(ImgStack_Reg_Input)-1)
    cuboid(4) = (width(ImgStack_Reg_Input)-1)-cuboid(1);
end

if (cuboid(2)+cuboid(5)) > (height(ImgStack_Reg_Input)-1)
    cuboid(5) = (height(ImgStack_Reg_Input)-1)-cuboid(2);
end

if (cuboid(3)+cuboid(6)) > (size(ImgStack_Reg_Input,3))
    cuboid(6) = (size(ImgStack_Reg_Input,3))-cuboid(3);
end

% = Crop
ImgStack_RegCrop_Output = imcrop3(ImgStack_Reg_Input, cuboid);


