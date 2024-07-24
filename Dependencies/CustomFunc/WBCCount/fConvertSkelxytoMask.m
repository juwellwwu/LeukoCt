function [Img_SkelPxFillMask_Out,Img_SkelPxFillMask_Area_Out] = fConvertSkelxytoMask(ImgStack_In,SkelPx_LinIdx_In,MaskAreaOnlyOption)

[SkelPx_x,SkelPx_y] = ind2sub([height(ImgStack_In) width(ImgStack_In)],SkelPx_LinIdx_In);
SkelPx_x = reshape(SkelPx_x,numel(SkelPx_x),1);
SkelPx_y = reshape(SkelPx_y,numel(SkelPx_y),1);

[SkelPx_xy_Boundary,Img_SkelPxFillMask_Area_Out] = boundary(SkelPx_x,SkelPx_y,0.75); 

% = Create Mask 
close all; 
if MaskAreaOnlyOption == 'n'
    Img_SkelPxFillMask_Out = zeros([height(ImgStack_In) width(ImgStack_In)],'logical');
    figure('Visible','on');
    imshow(Img_SkelPxFillMask_Out); 
    ROI_SkelPxFill = ... 
        images.roi.Polygon(gca,'Position',horzcat(SkelPx_y(SkelPx_xy_Boundary), SkelPx_x(SkelPx_xy_Boundary)));
    Img_SkelPxFillMask_Out = createMask(ROI_SkelPxFill);
else
    Img_SkelPxFillMask_Out = [];
end

close all;
