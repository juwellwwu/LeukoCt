function [ImgStack_Out,Img_Focus_Idx,Img_Focus_Idx_Illus] = ...
    fCreateImgStackFocusQualityMap(ImgStack_In,ImgStack_TissueMask_In,XY_PxLength_In)

CellDiameter_px_In = 5/XY_PxLength_In;

% = Ensure all slices of ImageStack has identical intensity, based on tissue pixels
ImgStack_ZIntensityScaleFactor = permute(sum(sum(ImgStack_In.*single(~ImgStack_TissueMask_In),1),2),[3,1,2]); 
ImgStack_ZIntensityScaleFactor = ImgStack_ZIntensityScaleFactor./min(ImgStack_ZIntensityScaleFactor,[],'all');

% Check
% figure;
% plot(ImgStack_ZIntensityScaleFactor);

ImgStack_Out = ...
    ImgStack_In./repmat(permute(ImgStack_ZIntensityScaleFactor,[2,3,1]),height(ImgStack_In),width(ImgStack_In));

% = Determine locations of best focus
% i) Z range (max-min)

ImgStack_Focus_Sm = imopen(ImgStack_Out,strel('disk',floor(0.5*CellDiameter_px_In))); 
ImgStack_Focus_Sm = ImgStack_Focus_Sm./max(ImgStack_Focus_Sm,[],'all');
ImgStack_Focus_Sm = ImgStack_Focus_Sm.*ImgStack_TissueMask_In;

ImgStack_Zmin = sort(ImgStack_Focus_Sm,3,"ascend"); 
ImgStack_Zmin = mean(ImgStack_Zmin(:,:,1:round(0.01*size(ImgStack_In,3))),3);
ImgStack_Zmin(~logical(ImgStack_TissueMask_In)) = NaN;

ImgStack_Zmax = sort(ImgStack_Focus_Sm,3,"descend");
ImgStack_Zmax = mean(ImgStack_Zmax(:,:,1:10),3);
ImgStack_Zmax(~logical(ImgStack_TissueMask_In)) = NaN;

ImgStack_Zrange = ImgStack_Zmax-ImgStack_Zmin;

% ii) Z projected standard deviation
ImgStack_Zstd = std(ImgStack_Focus_Sm,0,3);
ImgStack_Zstd(~logical(ImgStack_TissueMask_In)) = NaN;

% iii) Z projected variance
ImgStack_Zvar = var(ImgStack_Focus_Sm,0,3);
ImgStack_Zvar(~logical(ImgStack_TissueMask_In)) = NaN;

% clearvars ImgStack_Zmean;

% iii) Index ImgStack_Zrange_Idx and ImgStack_Zstd_Idx such that there are 10
% levels from high to low contrast
Zrange_Edge =  linspace(0.2,0.5,8); 
Img_Zrange_Idx = discretize(ImgStack_Zrange,Zrange_Edge)+1;
Img_Zrange_Idx(ImgStack_Zrange<0.2) = 1; 
Img_Zrange_Idx(ImgStack_Zrange>=0.5) = 10; 
Img_Zrange_Idx = single(Img_Zrange_Idx);

Zstd_Edge =  linspace(0.00,max(ImgStack_Zstd,[],'all','omitmissing'),10);
Img_Zstd_Idx = discretize(ImgStack_Zstd,Zstd_Edge);
Img_Zstd_Idx = single(Img_Zstd_Idx);

Zvar_Edge =  linspace(0.00,max(ImgStack_Zvar,[],'all','omitmissing'),10);
Img_Zvar_Idx = discretize(ImgStack_Zvar,Zvar_Edge);
Img_Zvar_Idx = single(Img_Zvar_Idx);

Img_Focus_Idx = 2/3*Img_Zrange_Idx+1/3*Img_Zstd_Idx; 

% Check: Prepare ImgStack_Focus_Idx_Illus
Img_Focus_Idx_Illus = imfuse(Img_Zrange_Idx./10,Img_Zstd_Idx./10,'ColorChannels',[1 2 0],'Scaling','none');  % pre 20230916
Img_Focus_Idx_Illus = horzcat(...
    255.*repmat(Img_Zrange_Idx./10,1,1,3),...
    255.*repmat(Img_Zstd_Idx./10,1,1,3),...
    Img_Focus_Idx_Illus,...
    255.*repmat(Img_Focus_Idx./10,1,1,3));

