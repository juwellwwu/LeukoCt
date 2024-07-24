function Img_LineROI = fCreateLineROIImg(LineROI_xlim,LineROI_ylim,ImgWidth,ImgHeight)

% Interpolate y values of Line ROI from xlim
LineROI_xlim_Lo = ceil(min(LineROI_xlim));
if LineROI_xlim_Lo <1
    LineROI_xlim_Lo = 1;
end
LineROI_xlim_Hi = floor(max(LineROI_xlim));
if LineROI_xlim_Hi > ImgWidth
    LineROI_xlim_Hi = ImgWidth;
end

LineROI_x_1 = LineROI_xlim_Lo:0.01:LineROI_xlim_Hi;
if numel(LineROI_x_1)>1
    LineROI_y_1 = interp1(LineROI_xlim',LineROI_ylim',LineROI_x_1,'linear');
    LineROI_y_1 = round(LineROI_y_1);
    LineROI_x_1 = round(LineROI_x_1);

    LineROI_x_1(LineROI_y_1<1) = [];
    LineROI_y_1(LineROI_y_1<1) = [];

    LineROI_x_1(LineROI_y_1>ImgHeight) = [];
    LineROI_y_1(LineROI_y_1>ImgHeight) = [];

    Img_LineROI_1 = zeros(ImgHeight,ImgWidth,'logical');
    Img_LineROI_1(sub2ind([height(Img_LineROI_1) width(Img_LineROI_1)],LineROI_y_1,LineROI_x_1)) = 1;
else
    Img_LineROI_1 = zeros(ImgHeight,ImgWidth,'logical');
end

% Interpolate x values of Line ROI from xlim
LineROI_ylim_Lo = ceil(min(LineROI_ylim));
if LineROI_ylim_Lo <1
    LineROI_ylim_Lo = 1;
end
LineROI_ylim_Hi = floor(max(LineROI_ylim));
if LineROI_ylim_Hi > ImgHeight
    LineROI_ylim_Hi = ImgHeight;
end

LineROI_y_2 = LineROI_ylim_Lo:0.01:LineROI_ylim_Hi;
if numel(LineROI_y_2)>1
    LineROI_x_2 = interp1(LineROI_ylim',LineROI_xlim',LineROI_y_2,'linear');
    LineROI_x_2 = round(LineROI_x_2);
    LineROI_y_2 = round(LineROI_y_2);

    LineROI_y_2(LineROI_x_2<1) = [];
    LineROI_x_2(LineROI_x_2<1) = [];

    LineROI_y_2(LineROI_x_2>ImgWidth) = [];
    LineROI_x_2(LineROI_x_2>ImgWidth) = [];

    Img_LineROI_2 = zeros(ImgHeight,ImgWidth,'logical');
    Img_LineROI_2(sub2ind([height(Img_LineROI_2) width(Img_LineROI_2)],LineROI_y_2',LineROI_x_2')) = 1;

else
    Img_LineROI_2 = zeros(ImgHeight,ImgWidth,'logical');
end

% Combine the two Img_LineROI, apply Img_AllSkel_Mask
Img_LineROI = max(Img_LineROI_1, Img_LineROI_2);
Img_LineROI = bwmorph(Img_LineROI,'skeleton');
