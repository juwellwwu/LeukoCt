function [Img_ThreshPerim_pCell,rp_ThreshPerim_pCell] = ...
    fCreateThreshPerim(Img_Thresh_Input,Img_ROI_Mask_Input)

Img_Temp_Perim = bwperim(Img_Thresh_Input);

% Remove edge pixels
Img_Temp_Perim(:,1:2) = 0;
Img_Temp_Perim(:,width(Img_Temp_Perim)-1:width(Img_Temp_Perim)) = 0;
Img_Temp_Perim(1:2,:) = 0;
Img_Temp_Perim(height(Img_Temp_Perim)-1:height(Img_Temp_Perim),:) = 0;

% Break perimeter w/ "spur' method
Img_Temp_Perim_Endpt = logical(Img_Temp_Perim.*~Img_ROI_Mask_Input);
rp_Endpt_Temp = regionprops(Img_Temp_Perim_Endpt,'Area','PixelIdxList');

Img_Temp_Perim_Endpt_Mask = zeros(height(Img_Temp_Perim_Endpt),width(Img_Temp_Perim_Endpt),height(rp_Endpt_Temp),'logical');
for epCt = 1:height(rp_Endpt_Temp) 
    Img_Clean = zeros(height(Img_Temp_Perim_Endpt),width(Img_Temp_Perim_Endpt),'logical');
    Img_Clean(rp_Endpt_Temp(epCt,1).PixelIdxList) = 1;
    if rp_Endpt_Temp(epCt,1).Area > 5
        Img_Clean = bwmorph(Img_Clean,'spur',round(0.5*rp_Endpt_Temp(epCt,1).Area-2));
    elseif rp_Endpt_Temp(epCt,1).Area < 3
        Img_Clean = imdilate(Img_Clean,strel('square',3));
    end
    Img_Temp_Perim_Endpt_Mask(:,:,epCt) = Img_Clean;
end

Img_Temp_Perim_Endpt_Mask = max(Img_Temp_Perim_Endpt_Mask,[],3);
Img_Temp_Perim = Img_Temp_Perim.*~Img_Temp_Perim_Endpt_Mask;

% Extract 2 longest lines for analysis
Img_Temp_Perim = ExtractNLargestBlobs(Img_Temp_Perim,2);

rp_Temp = regionprops(Img_Temp_Perim,'Area','PixelList','PixelIdxList');

Img_ThreshPerim_pCell = cell(2,1); 
rp_ThreshPerim_pCell = cell(2,1); 

% 1st edge
Img_ThreshPerim_pCell{1,1} = Img_Temp_Perim;
Img_ThreshPerim_pCell{1,1}(rp_Temp(1).PixelIdxList) = 0;
Img_ThreshPerim_pCell{1,1} = bwskel(Img_ThreshPerim_pCell{1,1}); 
Img_ThreshPerim_pCell{1,1} = Img_ThreshPerim_pCell{1,1}.*Img_ROI_Mask_Input;
Img_ThreshPerim_pCell{1,1}(1,:) = 0; 
Img_ThreshPerim_pCell{1,1}(end,:) = 0;
Img_ThreshPerim_pCell{1,1}(:,1) = 0;
Img_ThreshPerim_pCell{1,1}(:,end) = 0;
Img_ThreshPerim_pCell{1,1} = ExtractNLargestBlobs(Img_ThreshPerim_pCell{1,1},1); 
rp_ThreshPerim_pCell{1,1} = regionprops(Img_ThreshPerim_pCell{1,1},'Area','PixelList','PixelIdxList');

% 2nd edge
Img_ThreshPerim_pCell{2,1} = Img_Temp_Perim;
Img_ThreshPerim_pCell{2,1}(rp_Temp(2).PixelIdxList) = 0;
Img_ThreshPerim_pCell{2,1} = bwskel(Img_ThreshPerim_pCell{2,1}); 
Img_ThreshPerim_pCell{2,1} = Img_ThreshPerim_pCell{2,1}.*Img_ROI_Mask_Input;
Img_ThreshPerim_pCell{2,1}(1,:) = 0; 
Img_ThreshPerim_pCell{2,1}(end,:) = 0;
Img_ThreshPerim_pCell{2,1}(:,1) = 0;
Img_ThreshPerim_pCell{2,1}(:,end) = 0;
Img_ThreshPerim_pCell{2,1} = ExtractNLargestBlobs(Img_ThreshPerim_pCell{2,1},1); 
rp_ThreshPerim_pCell{2,1} = regionprops(Img_ThreshPerim_pCell{2,1},'Area','PixelList','PixelIdxList');