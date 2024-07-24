function [ImgStack_Crop_Out,...
    Img_Thresh_ROI_Crop_sCell_Out,rp_Thresh_ROI_Crop_sCell_Out,...
    Img_Skel_ROI_Crop_sCell_Out,rp_Skel_ROI_Crop_sCell_Out,...
    Img_ROI_Crop_Out,Img_Zstdev_Crop_Out,...
    PolygonROI_startXY_Out,ImgStack_ROICrop_cuboid,Img_ROI_Colocal_illus, ...
    ImgStack_In_Sz] = ...
    fPrepROIVesselImgStackParam(ImgStack_In,Img_ROI_In,Img_Zstdev_In,Img_Thresh_PreROI_sCell_In,Img_Skel_PreROI_sCell_In,PolygonROI_startXY_In,TimeParam_In)

rp_ROI = regionprops(Img_ROI_In,'Area','BoundingBox','PixelList','PixelIdxList'); 

%% Update Skeleton to only include segment in Img_ROI

Img_Skel_ROI_sCell = cell(2,1);
rp_Skel_ROI_sCell = cell(2,1);

for sCt = 1:2
    Img_Skel_ROI_sCell{sCt,1} = Img_Skel_PreROI_sCell_In{sCt,1}.*Img_ROI_In;
    rp_Skel_ROI_sCell{sCt,1} = regionprops(Img_Skel_ROI_sCell{sCt,1},'Area','PixelList','PixelIdxList'); 
end

% figure; 
% imshow(Img_ROI);

close all hidden;

clearvars Img_ROIDraw;
clearvars PolygonROI;
clearvars FigHandle;


%% Crop ImgStack to ROI area

ImgStack_ROICrop_cuboid = ... 
    [floor(rp_ROI.BoundingBox(1)) floor(rp_ROI.BoundingBox(2)) 1 (rp_ROI.BoundingBox(3)) (rp_ROI.BoundingBox(4)) size(ImgStack_In,3)-1];

if ImgStack_ROICrop_cuboid(1) < 1
    ImgStack_ROICrop_cuboid(1) = 1;
end

if ImgStack_ROICrop_cuboid(2) < 1
    ImgStack_ROICrop_cuboid(2) = 1;
end

while (round(ImgStack_ROICrop_cuboid(1))+round(ImgStack_ROICrop_cuboid(4))) > (width(ImgStack_In)-1)
    ImgStack_ROICrop_cuboid(4) = ImgStack_ROICrop_cuboid(4)-1;
end

while (round(ImgStack_ROICrop_cuboid(2))+round(ImgStack_ROICrop_cuboid(5))) > (height(ImgStack_In)-1)
    ImgStack_ROICrop_cuboid(5) = ImgStack_ROICrop_cuboid(5)-1;
end

% = Crop ImgStack_In
ImgStack_Crop_Out = imcrop3(ImgStack_In, ImgStack_ROICrop_cuboid);
ImgStack_Crop_Out = ImgStack_Crop_Out(:,:,1:floor(size(ImgStack_Crop_Out,3)/TimeParam_In.tSeg_FrLength).*TimeParam_In.tSeg_FrLength);

fprintf(strcat('\nHeight of ImgStack_Crop =',32,num2str(height(ImgStack_Crop_Out),'%0d'),'.\n'));
fprintf(strcat('Width of ImgStack_Crop =',32,num2str(width(ImgStack_Crop_Out),'%0d'),'.\n'));
fprintf(strcat('Depth of ImgStack_Crop =',32,num2str(size(ImgStack_Crop_Out,3),'%0d'),'.\n\n'));


%% Crop Img_Skel, Img_Thresh, Img_ROI, Img_ZStdev to ROI area

fprintf('Cropping Img_Skel, Img_Thresh, Img_ROI to ROI area...\n');

% = Crop Threshold & Skeleton Image
Img_Thresh_ROI_Crop_sCell_Out = cell(2,1);
Img_Skel_ROI_Crop_sCell_Out = cell(2,1);

rp_Thresh_ROI_Crop_sCell_Out = cell(2,1);
rp_Skel_ROI_Crop_sCell_Out = cell(2,1);

for sCt = 1:2
    Img_Thresh_ROI_Crop_sCell_Out{sCt,1} = imcrop(Img_Thresh_PreROI_sCell_In{sCt,1},ImgStack_ROICrop_cuboid([1,2,4,5]));
    rp_Thresh_ROI_Crop_sCell_Out{sCt,1} = regionprops(Img_Thresh_ROI_Crop_sCell_Out{sCt,1},'Area','PixelList','PixelIdxList'); 

    Img_Skel_ROI_Crop_sCell_Out{sCt,1} = imcrop(Img_Skel_PreROI_sCell_In{sCt,1},ImgStack_ROICrop_cuboid([1,2,4,5]));
    rp_Skel_ROI_Crop_sCell_Out{sCt,1} = regionprops(Img_Skel_ROI_Crop_sCell_Out{sCt,1},'Area','PixelList','PixelIdxList'); 
end

Img_ROI_Crop_Out = imcrop(Img_ROI_In, ImgStack_ROICrop_cuboid([1,2,4,5]));
rp_ROI_Crop_Out = regionprops(Img_ROI_Crop_Out,'Area','PixelList','PixelIdxList'); % Single blob from LineROI

Img_Zstdev_Crop_Out = imcrop(Img_Zstdev_In,ImgStack_ROICrop_cuboid([1,2,4,5]));


%% Update Img_Thresh_ROI_Crop & Img_Skel_ROI_Crop: Keep single blob or skeleton location matched to vessel ROI

% == Find which blob in Img_Thresh_ROI_Crop & Img_Skel_ROI_Crop is vessel ROI
% = Thresh
Thresh_Vessel_BlobIdx_sCell = cell(2,1);
for sCt=1:2
    min_pDist_Thresh =zeros(size(rp_Thresh_ROI_Crop_sCell_Out{1,1},1),1);
    for bCt=1:size(rp_Thresh_ROI_Crop_sCell_Out{sCt,1},1) % blob count
        min_pDist_Thresh(bCt,1) = ...
            mean(pdist2(rp_ROI_Crop_Out.PixelList,rp_Thresh_ROI_Crop_sCell_Out{sCt,1}(bCt).PixelList,'euclidean','Smallest',1));
    end
    [~,Thresh_Vessel_BlobIdx_sCell{sCt,1}] = min(min_pDist_Thresh); 
end

% = Skel
Skel_Vessel_BlobIdx_sCell = cell(2,1);
for sCt = 1:2
    min_pDist_Skel =zeros(size(rp_Skel_ROI_Crop_sCell_Out{sCt,1},1),1);
    for bCt=1:size(rp_Skel_ROI_Crop_sCell_Out{sCt,1},1) 
        min_pDist_Skel(bCt,1) = ...
            mean(pdist2(rp_ROI_Crop_Out.PixelList,rp_Skel_ROI_Crop_sCell_Out{sCt,1}(bCt).PixelList,'euclidean','Smallest',1));
    end
    [~,Skel_Vessel_BlobIdx_sCell{sCt,1}] = min(min_pDist_Skel); 
end


% == Update Img_Thresh_ROI_Crop & Img_Skel_ROI_Crop
for sCt = 1:2
    Img_Thresh_ROI_Crop_sCell_Out{sCt,1} = zeros(size(Img_Thresh_ROI_Crop_sCell_Out{sCt,1}),'logical');
    Img_Thresh_ROI_Crop_sCell_Out{sCt,1}(rp_Thresh_ROI_Crop_sCell_Out{sCt,1}(Thresh_Vessel_BlobIdx_sCell{sCt,1}).PixelIdxList) = 1;
    rp_Thresh_ROI_Crop_sCell_Out{sCt,1} = regionprops(Img_Thresh_ROI_Crop_sCell_Out{sCt,1},'Area','PixelList','PixelIdxList'); % Update

    Img_Skel_ROI_Crop_sCell_Out{sCt,1} = zeros(size(Img_Skel_ROI_Crop_sCell_Out{sCt,1}),'logical');
    Img_Skel_ROI_Crop_sCell_Out{sCt,1}(rp_Skel_ROI_Crop_sCell_Out{sCt,1}(Skel_Vessel_BlobIdx_sCell{sCt,1}).PixelIdxList) = 1;
    rp_Skel_ROI_Crop_sCell_Out{sCt,1} = regionprops(Img_Skel_ROI_Crop_sCell_Out{sCt,1},'Area','PixelList','PixelIdxList'); 
end

clearvars Thresh_Vessel_BlobIdx_sCell min_pDist_Thresh sCt bCt;
clearvars Skel_Vessel_BlobIdx_sCell min_pDist_Skel sCt bCt;


%%  Move PolygonROI_startXY based on ROI 

PolygonROI_startXY_Out = PolygonROI_startXY_In - [rp_ROI.BoundingBox(1) rp_ROI.BoundingBox(2)];


%%  Generate image showing location of vessel ROI, Img_Thresh_ROI_Crop

figure; 
Img_ROI_Colocal_illus = horzcat(255.*repmat(min(ImgStack_Crop_Out,[],3),1,1,3),...
    255.*repmat(Img_Zstdev_Crop_Out,1,1,3), ...
    imfuse(Img_Thresh_ROI_Crop_sCell_Out{1,1}.*0.50.*Img_ROI_Crop_Out,repmat(min(ImgStack_Crop_Out,[],3),1,1,3),'ColorChannels',[1 2 0],'Scaling','none'),...
    imfuse(Img_Thresh_ROI_Crop_sCell_Out{2,1}.*0.50.*Img_ROI_Crop_Out,Img_Zstdev_Crop_Out,'ColorChannels',[1 2 0],'Scaling','none'),...
    imfuse(min(ImgStack_Crop_Out,[],3),Img_Zstdev_Crop_Out,'ColorChannels',[1 2 0],'Scaling','none'),...
    imfuse(Img_Thresh_ROI_Crop_sCell_Out{1,1}.*0.50.*Img_ROI_Crop_Out,Img_Thresh_ROI_Crop_sCell_Out{2,1}.*0.50.*Img_ROI_Crop_Out,'ColorChannels',[1 2 0],'Scaling','none'));
imshow(Img_ROI_Colocal_illus,'InitialMagnification',200);

% close all;

clearvars Thresh*_Crop_BlobCt Skel*_Crop_BlobCt;
clearvars Thresh*_Vessel_BlobIdx Skel*_Vessel_BlobIdx;
clearvars ROIThresh*_MinDist ROISkel*_MinDist;
clearvars Img_Temp;


%% Save size of ImgStack_In
% Needed for creating Demo Illustration Image

ImgStack_In_Sz = size(ImgStack_In);

