function [Img_Zstdev_Out,Img_Thresh_PreROI_PreClean_sCell_Out,Img_Thresh_PreROI_sCell_Out,Img_Skel_PreROI_sCell_Out] = ...
    fCreatePreROIDrawVesselMask(ImgStack_In,XY_PxLength_In,SaveFilePath_In)

VesselLumenOnlyOption = 'n';
ExtraVesselAreaOption = 'y'; % Default: 'n'; use 'y' if masks are missing chunks 
PlotOption = 'n';

Img_Thresh_PreROI_PreClean_sCell_Out = ...
   repmat({zeros(height(ImgStack_In),width(ImgStack_In),2)},2,1); 

[Img_Zstdev_Out,Img_Thresh_PreROI_PreClean_sCell_Out{1,1},Img_Thresh_PreROI_PreClean_sCell_Out{2,1}] = ...
    fCreateRegImgZthresh(ImgStack_In,XY_PxLength_In,SaveFilePath_In,VesselLumenOnlyOption,ExtraVesselAreaOption,PlotOption);

% = Check
% figure; 
% Img_Temp = imfuse(Img_Thresh_PreClean_sCell_Out{1,1}.*255,Img_Thresh_PreClean_sCell_Out{2,1}.*255,'Scaling','none','ColorChannels',[1 2 0]);
% imshow(Img_Temp,'InitialMagnification',200);

clearvars k Img_Temp;


%%
% == Crude mask cleanup, step 1: user Identifies vessel to keep 
% Only region in approx mask defined as foreground (blue) is kept

GraphCut_ROI_Stack = zeros(size(Img_Thresh_PreROI_PreClean_sCell_Out{1,1},1),size(Img_Thresh_PreROI_PreClean_sCell_Out{1,1},2),3);

FigHandle = figure;
set(FigHandle, 'MenuBar', 'none');
set(FigHandle, 'ToolBar', 'none');
imshow(Img_Thresh_PreROI_PreClean_sCell_Out{1,1}.*Img_Zstdev_Out,'InitialMagnification',300);

% Interactive: draw foreground and background lines
% [x,y] positions of points on line in ROI.Position
fprintf('**USER INPUT**: draw foreground poly-line; poly-line should be within vessel to be kept...\n');
fprintf('(Double-click to end draw)...\n');
fgd_ROI = drawpolyline;
GraphCut_ROI_Stack(:,:,1) = createMask(fgd_ROI);

% All 0 px count as background
GraphCut_ROI_Stack(:,:,2) = ~Img_Thresh_PreROI_PreClean_sCell_Out{1,1}; 

% Keep only blobs containing fgd_ROI.Position
Img_Thresh = Img_Thresh_PreROI_PreClean_sCell_Out{1,1};

rp_Thresh_Temp = regionprops(Img_Thresh,'PixelIdxList');
for bCt = 1:numel(rp_Thresh_Temp)
    isMember = ... 
        logical(round(mean(ismember(sub2ind(size(Img_Thresh_PreROI_PreClean_sCell_Out{1,1}),round(fgd_ROI.Position(:,2)),round(fgd_ROI.Position(:,1))),rp_Thresh_Temp(bCt).PixelIdxList))));
    if isMember<1 
        Img_Thresh(rp_Thresh_Temp(bCt).PixelIdxList) = 0;
    end
end

% Elect largest blob
% Function ExtractNLargestBlobs() by user "Image Analyst" in MATLAB community: 
% https://www.mathworks.com/matlabcentral/answers/68696-how-can-i-extract-the-largest-blob-in-a-binary-image#answer_79991
Img_Thresh = ExtractNLargestBlobs(Img_Thresh,1);

clearvars fgd_ROI bgd_ROI break_ROI;
clearvars rp_Thresh_Temp isMember bCt Img_Thresh_Colocal;


%%
% == Crude mask cleanup, step 2: user applies line eraser
% Allow user to trim Img_Thresh, clean up vessel boundaries

FigHandle = figure;
set(FigHandle, 'MenuBar', 'none');
set(FigHandle, 'ToolBar', 'none');
Img_Thresh_Colocal = imfuse(0.33*Img_Thresh,Img_Zstdev_Out,'Scaling','none','ColorChannels',[1 2 0]);
imshow(Img_Thresh_Colocal,'InitialMagnification',300);

fprintf('**USER INPUT**: draw line eraser...\n');
fprintf('*Only largest blob will be kept. Double-click to end draw ...\n');
break_ROI = drawpolyline('Color',[0.4660 0.6740 0.1880]);
GraphCut_ROI_Stack(:,:,3) = imdilate(createMask(break_ROI),strel('disk',3));

clearvars Img_Thresh_Colocal;


%%
% == Finalize cleaned up Img_Thresh (no user input)

Img_Thresh_PreROI_sCell_Out = ... 
    repmat({zeros(height(Img_Thresh),width(Img_Thresh),2,'logical')},2,1);

% i) Vessel Lumen Mask
Img_Thresh = and(Img_Thresh,~GraphCut_ROI_Stack(:,:,3));
Img_Thresh = ExtractNLargestBlobs(Img_Thresh,1); 

sgolay_FrLength = 2*round(0.5*12/XY_PxLength_In)+1; 
[~,Img_Thresh] = fSmoothImgBWBoundary_SG(Img_Thresh,2,sgolay_FrLength,XY_PxLength_In);  
Img_Thresh = ExtractNLargestBlobs(Img_Thresh,1);
Img_Thresh_PreROI_sCell_Out{1,1} = Img_Thresh;

% ii) Flow profile Mask
Img_Thresh_PreROI_sCell_Out{2,1} = Img_Thresh_PreROI_PreClean_sCell_Out{2,1};
[~,Img_Thresh_PreROI_sCell_Out{2,1}] = fSmoothImgBWBoundary_SG(Img_Thresh_PreROI_sCell_Out{2,1},2,sgolay_FrLength,XY_PxLength_In);  
Img_Thresh_PreROI_sCell_Out{2,1} = logical(Img_Thresh_PreROI_sCell_Out{2,1}.*Img_Thresh_PreROI_sCell_Out{1,1}); 
Img_Thresh_PreROI_sCell_Out{2,1} = ExtractNLargestBlobs(Img_Thresh_PreROI_sCell_Out{2,1},1);

% figure; 
% imshowpair(Img_Zstdev,Img_Thresh);


%%
% == Calculate skeletons for pre-ROI Vessel Masks
Img_Skel_PreROI_sCell_Out = ... 
    repmat({zeros(size(Img_Thresh_PreROI_sCell_Out),'logical')},2,1);

for sCt = 1:2
    Img_Skel_PreROI_sCell_Out{sCt,1} = ...
        fCreateImgSkel(Img_Thresh_PreROI_sCell_Out{sCt,1},ceil(1.2*max(bwdist(bwperim(Img_Thresh_PreROI_sCell_Out{sCt,1})).*Img_Thresh_PreROI_sCell_Out{sCt,1},[],'all')));
end


% figure;
% imshowpair(Img_Skel,Img_Zstdev,'falseColor');
