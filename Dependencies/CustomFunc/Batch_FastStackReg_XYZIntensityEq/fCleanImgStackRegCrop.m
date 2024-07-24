function [ImgStack_Clean_Output] = fCleanImgStackRegCrop(ImgStack_Input)

HorzLine=sum(diff(range(ImgStack_Input,3),1,1),2); 
VertLine=sum(diff(range(ImgStack_Input,3),1,2),1); 

HorzLine_Outlier = isoutlier(HorzLine,'median',"ThresholdFactor",15);
HorzLine_Outlier = find(HorzLine_Outlier);
HorzLine_Outlier_Top = max(HorzLine_Outlier(HorzLine_Outlier<0.25*height(ImgStack_Input)),[],'all');
HorzLine_Outlier_Bttm = min(HorzLine_Outlier(HorzLine_Outlier>0.75*height(ImgStack_Input)),[],'all')-1; 

VertLine_Outlier = isoutlier(VertLine,'median',"ThresholdFactor",15);
VertLine_Outlier = find(VertLine_Outlier);
VertLine_Outlier_Left = max(VertLine_Outlier(VertLine_Outlier<0.25*width(ImgStack_Input)),[],'all');
VertLine_Outlier_Right = min(VertLine_Outlier(VertLine_Outlier>0.75*width(ImgStack_Input)),[],'all'); 

% Create Mask describing the borders to remove
Mask = zeros(height(ImgStack_Input),width(ImgStack_Input),'logical');

if ~isempty(HorzLine_Outlier_Top)
    Mask(1:HorzLine_Outlier_Top,:)  = 1;
end

if ~isempty(HorzLine_Outlier_Bttm)
    Mask(HorzLine_Outlier_Bttm:end,:)  = 1;
end

if ~isempty(VertLine_Outlier_Left)
    Mask(:,1:VertLine_Outlier_Left)  = 1;
end

if ~isempty(VertLine_Outlier_Right)
    Mask(:,VertLine_Outlier_Right:end)  = 1;
end

Mask = imcomplement(Mask);

% Check
% % % figure;
% % % imshow(Mask);

rp_Mask = regionprops(Mask,'BoundingBox');

ImgStack_Clean_Output = ...
    imcrop3(ImgStack_Input,[ceil(rp_Mask(1).BoundingBox(1)) ceil(rp_Mask(1).BoundingBox(2)) 1 floor(rp_Mask(1).BoundingBox(3))-1 floor(rp_Mask(1).BoundingBox(4))-1 size(ImgStack_Input,3)-1]);


