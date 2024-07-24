function [Img_Zstdev_Output,Img_Mask_VesselLumen,Img_Mask_Flow] = ...
    fCreateRegImgZthresh(ImgStack_Input,XY_PxLength_Input,SaveFilePath_Input,VesselLumenOnlyOption,ExtraVesselAreaOption,PlotOption)


CellDiameter_px_Input = 5/XY_PxLength_Input;


%% Locate pixels w/out significant plane-to-plane intensity changes

ImgStack_Zmin_Mask = imbinarize(min(ImgStack_Input,[],3)); 

absDiffLimit = (0.050:0.005:0.15)'; % default: (0.05:0.005:0.15)'
ImgStack_absDiff = ... 
    ImgStack_Input.*imclose(~ImgStack_Zmin_Mask,strel('disk',round(1.00*CellDiameter_px_Input))); 
ImgStack_absDiff = abs(diff(ImgStack_absDiff,1,3));

absDiff_NaNIdx_Ct = zeros(numel(absDiffLimit),1);
for lCt = 1:numel(absDiffLimit) 
    absDiff_NaNIdx_Ct(lCt,1) = sum(ImgStack_absDiff<absDiffLimit(lCt,1),'all');
end

absDiff_NaNIdx_Ct_Sm = absDiff_NaNIdx_Ct;
for smCt = 1:3 
    absDiff_NaNIdx_Ct_SmObj = fit((1:1:numel(absDiffLimit))',double(absDiff_NaNIdx_Ct_Sm),'smoothingspline');
    absDiff_NaNIdx_Ct_Sm = feval(absDiff_NaNIdx_Ct_SmObj,(1:1:numel(absDiffLimit))');
end

absDiffLimitIdx = knee_pt(absDiff_NaNIdx_Ct_Sm,(1:1:numel(absDiffLimit))',1);

if PlotOption == 'y'
    FigHandle0 = figure('Position',[10 700 300 200],'Visible','off'); 
    tlo = tiledlayout(1,1,'TileSpacing','Compact','Padding','tight');
    ax = nexttile(tlo);
    hold(ax, 'on');
    plot(absDiff_NaNIdx_Ct,'LineWidth',1);
    if exist('absDiff_NaNIdx_Ct_Sm')
        plot(absDiff_NaNIdx_Ct_Sm,':','LineWidth',2);
    end
    xline(absDiffLimitIdx);
    xlabel('absDiffLimit','interpreter','none');
    ylabel('absDiff_NaNIdx_Ct','interpreter','none');
    title(strcat('Chosen absDiffLimitIdx =',32,num2str(absDiffLimitIdx),'; absDiffLimit =',32,num2str(absDiffLimit(absDiffLimitIdx))),'Interpreter','none');
    fontsize(gcf,6,"points");
    saveas(FigHandle0,strcat(SaveFilePath_Input,'absDiff_NaNIdx_Ct_Plot.tif'));
end

ImgStack_absDiff = abs(diff(ImgStack_Input,1,3)); 
absDiffLimit = absDiffLimit(absDiffLimitIdx);
ImgStack_absDiff_NaNIdx = ImgStack_absDiff<absDiffLimit; 
ImgStack_absDiff = cat(3,zeros(height(ImgStack_absDiff),width(ImgStack_absDiff),'single'),ImgStack_absDiff);


%% Moving Z standard deviation

Img_Zstdev_Output = ImgStack_Input;
Img_Zstdev_Output(ImgStack_absDiff_NaNIdx) = NaN; 
Img_Zstdev_Output = std(ImgStack_Input,0,3,'omitnan'); 

Img_Zstdev_MinIntensity = min(Img_Zstdev_Output,[],'all');
Img_Zstdev_MaxIntensity = max(Img_Zstdev_Output,[],'all');
Img_Zstdev_Output = (Img_Zstdev_Output-Img_Zstdev_MinIntensity).*(1-0)./(Img_Zstdev_MaxIntensity-Img_Zstdev_MinIntensity);


%% Define approx tissue mask as fulfilling both criteria defined by ImgStack_ZMeanBin and ImgStack_Zmin_Mask

Img_Tissue = logical(min(~imbinarize(Img_Zstdev_Output),ImgStack_Zmin_Mask));
Img_NonFlow = logical(~imbinarize(Img_Zstdev_Output));

% % % if PlotOption == 'y'
% % %     figure;
% % %     Img = imshow(horzcat(...
% % %         255.*repmat(~Img_Tissue,1,1,3),...
% % %         255.*repmat(~Img_NonFlow,1,1,3),...
% % %         imfuse(~Img_Tissue,~Img_NonFlow,'ColorChannels',[1 2 0],'Scaling','none')),...
% % %         'InitialMagnification',200);
% % %     Img = Img.CData;
% % %     imwrite(Img,strcat(SaveFilePath_Input,'ImgTissue_ImgNonFlow.tif'));
% % % end

close all;
clearvars Img;


%% Create Vessel Lumen & Flow Mask based on absDiff 

ImgStack_absDiff = single(logical(ImgStack_absDiff));
ImgStack_absDiff(ImgStack_absDiff_NaNIdx) = NaN;
ImgStack_absDiffSum = sum(ImgStack_absDiff,3,'omitnan'); 


%% Create Vessel Lumen & Flow Mask based on absDiff 

Frac = (0.002:0.002:0.20); % default (0.002:0.002:0.20)
ImgStack_absDiffSum_Mask = zeros(height(ImgStack_absDiffSum),width(ImgStack_absDiffSum),numel(Frac),'logical'); 
Tissue_absDiffSum_ColocalPxFrac = zeros(numel(Frac),1,'single'); 
NonFlow_absDiffSum_ColocalPxFrac = zeros(numel(Frac),1,'single'); 

for fracCt = 1:numel(Frac)

    Img_Temp = ImgStack_absDiffSum_Mask(:,:,fracCt); 
    Img_Temp(ImgStack_absDiffSum < Frac(fracCt)*size(ImgStack_Input,3)) = 0; 
    Img_Temp(ImgStack_absDiffSum >= Frac(fracCt)*size(ImgStack_Input,3)) = 1; 

    ImgStack_absDiffSum_Mask(:,:,fracCt) = logical(Img_Temp);

    Tissue_absDiffSum_ColocalPxFrac(fracCt,1) = ... 
        sum(~Img_Temp.*Img_Tissue,'all')/sum(Img_Tissue,'all');
    NonFlow_absDiffSum_ColocalPxFrac(fracCt,1) = ... 
        sum(~Img_Temp.*Img_NonFlow,'all')/sum(Img_NonFlow,'all');

end

clearvars Img_Temp fracCt;


%% Create Vessel Lumen & Flow Mask based on absDiff 

Tissue_absDiffSum_ColocalPxFrac_Sm = Tissue_absDiffSum_ColocalPxFrac;
NonFlow_absDiffSum_ColocalPxFrac_Sm = NonFlow_absDiffSum_ColocalPxFrac;

for smCt = 1:3 

    Tissue_absDiffSum_ColocalPxFrac_SmObj = fit((1:1:numel(Frac))',double(Tissue_absDiffSum_ColocalPxFrac_Sm),'smoothingspline');
    Tissue_absDiffSum_ColocalPxFrac_Sm = feval(Tissue_absDiffSum_ColocalPxFrac_SmObj,(1:1:numel(Frac))');

    NonFlow_absDiffSum_ColocalPxFrac_SmObj = fit((1:1:numel(Frac))',double(NonFlow_absDiffSum_ColocalPxFrac_Sm),'smoothingspline');
    NonFlow_absDiffSum_ColocalPxFrac_Sm = feval(NonFlow_absDiffSum_ColocalPxFrac_SmObj,(1:1:numel(Frac))');

end

if ExtraVesselAreaOption == 'n' 
    VesselLumenMaskFracIdx1 = knee_pt(Tissue_absDiffSum_ColocalPxFrac_Sm,(1:1:numel(Frac))',1);
    [~,VesselLumenMaskFracIdx2] = knee_pt((1:1:numel(Frac))',Tissue_absDiffSum_ColocalPxFrac_Sm,1); 
    VesselLumenMaskFracIdx = round(0.5*(VesselLumenMaskFracIdx1+VesselLumenMaskFracIdx2));
else 
    VesselLumenMaskFracIdx = knee_pt(Tissue_absDiffSum_ColocalPxFrac_Sm,(1:1:numel(Frac))',1);
    VesselLumenMaskFracIdx1 = VesselLumenMaskFracIdx;
    VesselLumenMaskFracIdx2 = VesselLumenMaskFracIdx;
end

if VesselLumenOnlyOption == 'n' 
    if ExtraVesselAreaOption == 'n'
        FlowMaskFracIdx1 = knee_pt(NonFlow_absDiffSum_ColocalPxFrac_Sm,(1:1:numel(Frac))',1);
        [~,FlowMaskFracIdx2] = knee_pt((1:1:numel(Frac))',NonFlow_absDiffSum_ColocalPxFrac_Sm,1); 
    else
        FlowMaskFracIdx1 = knee_pt(NonFlow_absDiffSum_ColocalPxFrac_Sm,(1:1:numel(Frac))',1);
        FlowMaskFracIdx2 = FlowMaskFracIdx1;
    end

    FlowMaskFracIdx = ...
        find(movmean(diff(NonFlow_absDiffSum_ColocalPxFrac)./(NonFlow_absDiffSum_ColocalPxFrac(1:end-1))*100,3)<0.1)+1; 
    FlowMaskFracIdx = FlowMaskFracIdx(FlowMaskFracIdx>FlowMaskFracIdx2);
    FlowMaskFracIdx = min(FlowMaskFracIdx,[],'all');
    if isempty(FlowMaskFracIdx)
        FlowMaskFracIdx = numel(Frac);
    end
end


%% Figures, Vessel Lumen Mask based on absDiff 
% Do not switch order of this section with sections afterwards

if PlotOption == 'y'

    % i. Different Tissue Frac Colocalized with absDiff
    figure;
    Img1 = montage(ImgStack_absDiffSum_Mask.*0.6,'Size',[NaN,10]); 
    Img1 = Img1.CData;
    Img2 = montage(repmat(min(ImgStack_Input,[],3),1,1,numel(Frac)),'Size',[NaN,10]); 
    Img2 = Img2.CData;
    Img = imfuse(Img1,Img2,'ColorChannels',[1 2 0]);
    imwrite(Img,strcat(SaveFilePath_Input,'Tissue_absDiffSumMask_Colocal_Montage_ChosenFrac',num2str(VesselLumenMaskFracIdx),'.tif'));
    close all;

    % ii. Knee-point of Tissue_absDiffSum_ColocalPxFrac
    FigHandle1 = figure('Position',[10 400 300 200],'Visible','off'); 
    tlo = tiledlayout(1,1,'TileSpacing','Compact','Padding','tight');
    ax = nexttile(tlo);
    hold(ax, 'on');
    plot(100*Tissue_absDiffSum_ColocalPxFrac,'LineWidth',1,'Color',[0.5 0 0]);
    if exist('Tissue_absDiffSum_ColocalPxFrac_Sm')
        plot(100*Tissue_absDiffSum_ColocalPxFrac_Sm,':','LineWidth',2,'Color',[0.5 0 0]);
    end
    xline(VesselLumenMaskFracIdx,'LineWidth',0.5,'Color',[0.25 0 0]);
    if exist('VesselLumenMaskFracIdx1')
        xline(VesselLumenMaskFracIdx1,':','LineWidth',0.5,'Color',[0.25 0 0]);
    end
    if exist('VesselLumenMaskFracIdx2')
        xline(VesselLumenMaskFracIdx2,':','LineWidth',0.5,'Color',[0.25 0 0]);
    end
    xlabel('Frac Index','interpreter','none');
    ylabel({'Colocalized Pixel %:';'Tissue, Area excluded by absDiffSumMask'},'interpreter','none');
    title(strcat('Chosen Frac Index =',32,num2str(VesselLumenMaskFracIdx),'; Frac =',32,num2str(Frac(VesselLumenMaskFracIdx))),'interpreter','none');
    fontsize(gcf,6,"points");
    saveas(FigHandle1,strcat(SaveFilePath_Input,'Tissue_absDiffSum_ColocalPxFrac_Plot.tif'));

end

clearvars Img1 Img2 Img;


%% Figures, Flow Mask based on absDiff 

if VesselLumenOnlyOption == 'n' & PlotOption == 'y'

    % i. Different NonFlowArea Frac Colocalized with absDiff
    figure;
    Img1 = montage(ImgStack_absDiffSum_Mask.*0.6,'Size',[NaN,10]); 
    Img1 = Img1.CData;
    Img2 = montage(repmat(Img_Zstdev_Output,1,1,numel(Frac)),'Size',[NaN,10]); 
    Img2 = Img2.CData;
    Img = imfuse(Img1,Img2,'ColorChannels',[1 2 0]);
    imwrite(Img,strcat(SaveFilePath_Input,'NonFlow_absDiffSumMask_Colocal_Montage_ChosenFrac',num2str(FlowMaskFracIdx),'.tif'));
    close all;

    % ii. Knee-point of NonFlow_absDiffSum_ColocalPxFrac
    FigHandle1 = figure('Position',[10 400 300 200],'Visible','off'); 
    tlo = tiledlayout(1,1,'TileSpacing','Compact','Padding','tight');
    ax = nexttile(tlo);
    hold(ax, 'on');
    plot(100*NonFlow_absDiffSum_ColocalPxFrac,'LineWidth',1,'Color',[0 0 0.5]);
    if exist('NonFlow_absDiffSum_ColocalPxFrac_Sm')
        plot(100*NonFlow_absDiffSum_ColocalPxFrac_Sm,':','LineWidth',2,'Color',[0 0 0.5]);
    end
    xline(FlowMaskFracIdx,'LineWidth',0.5,'Color',[0 0 0.25]);
    if exist('FlowMaskFracIdx1')
        xline(FlowMaskFracIdx1,':','LineWidth',0.5,'Color',[0 0 0.25]);
    end
    if exist('FlowMaskFracIdx2')
        xline(FlowMaskFracIdx2,':','LineWidth',0.5,'Color',[0 0 0.25]);
    end
    xlabel('Frac Index','interpreter','none');
    ylabel({'Colocalized Pixel %:';'NonFlowArea, Area excluded by absDiffSumMask'},'interpreter','none');
    title(strcat('Chosen Frac Index =',32,num2str(FlowMaskFracIdx),'; Frac =',32,num2str(Frac(FlowMaskFracIdx))),'interpreter','none');
    fontsize(gcf,6,"points");
    saveas(FigHandle1,strcat(SaveFilePath_Input,'NonFlowArea_absDiffSum_ColocalPxFrac_Plot.tif'));

end

clearvars Img1 Img2 Img;


%% Create Mask based on absDiff to clean tissue background

% = Elect Vessel Lumen Mask using VesselLumenMaskFracIdx
if ~isnan(VesselLumenMaskFracIdx)
    Img_Mask_VesselLumen = ImgStack_absDiffSum_Mask(:,:,VesselLumenMaskFracIdx);
else
    fprintf('... Kneepoint of Tissue_absDiffSum_ColocalPxFrac cannot be found. Assume Frac = 0.05.\n');
    Img_Mask_VesselLumen = ImgStack_absDiffSum_Mask(:,:,Frac==0.05);    
end

% = Cleanup
Img_Mask_VesselLumen = bwmorph(Img_Mask_VesselLumen,'majority');

Area_Limit = (5/XY_PxLength_Input).^2.*5;
rp_Mask_VesselLumen = regionprops(Img_Mask_VesselLumen,'Area','PixelIdxList');
rp_Mask_VesselLumen_Area = cat(1,rp_Mask_VesselLumen.Area);
rp_Mask_VesselLumen_Area_FlagID = find(rp_Mask_VesselLumen_Area<Area_Limit);
Img_Mask_VesselLumen(cat(1,rp_Mask_VesselLumen(rp_Mask_VesselLumen_Area_FlagID).PixelIdxList)) = 0;


%% Create Mask based on absDiff to measure flow velocity

if VesselLumenOnlyOption == 'n' 
    % = Elect Mask using FlowVelocityMaskFracIdx
    Img_Mask_Flow = ImgStack_absDiffSum_Mask(:,:,FlowMaskFracIdx);

    % = Cleanup
    Img_Mask_Flow = bwmorph(Img_Mask_Flow,'majority');
end


%% Clean up Masks for Output

Img_Mask_VesselLumen = ... 
    logical(Img_Mask_VesselLumen.* ...
    imdilate(~Img_Tissue,strel('disk',round(1.0*CellDiameter_px_Input)))); 
rp_Mask_VesselLumen = regionprops(Img_Mask_VesselLumen,'PixelIdxList');

% Smooth boundary of thresholded blob area
[Img_Mask_VL_SmoothBoundary] = ...
    fSmoothImgBWBoundary_SG(Img_Mask_VesselLumen,2,2*round(0.5*1.0*CellDiameter_px_Input)+1,XY_PxLength_Input); % default: order=2; frLength = 9;
rp_Mask_VL_SmoothBoundary = regionprops(Img_Mask_VL_SmoothBoundary,'PixelIdxList');
rp_Mask_VL_SmoothBoundary_PixelIdxList = cat(1,rp_Mask_VL_SmoothBoundary.PixelIdxList);

% Remove blobs from mask skipped by fSmoothImgBWBoundary_SG()
FlagID = zeros(height(rp_Mask_VesselLumen),1,'logical');
for BlobCt = 1:height(rp_Mask_VesselLumen)
    FlagID(BlobCt) = ~logical(numel(find(ismember(rp_Mask_VL_SmoothBoundary_PixelIdxList,rp_Mask_VesselLumen(BlobCt).PixelIdxList))));
end
Img_Mask_VesselLumen(cat(1,rp_Mask_VesselLumen(FlagID).PixelIdxList)) = 0;

% Fill residual small holes in mask, if exist
Img_Mask_VesselLumen = ~bwareaopen(~Img_Mask_VesselLumen,round(CellDiameter_px_Input.^2));

% Fill region within smoothed boundary; update Img_Mask1
Img_Mask_VesselLumen_Skel = bwmorph(bwskel(Img_Mask_VesselLumen),'spur',3);
Img_Mask_VesselLumen_Skel = Img_Mask_VesselLumen_Skel.*~Img_Mask_VL_SmoothBoundary;
rp_Mask_VesselLumen_Skel = regionprops(Img_Mask_VesselLumen_Skel,'Area','PixelIdxList');
FlagID = cat(1,rp_Mask_VesselLumen_Skel.Area) < (CellDiameter_px_Input*4); % expect skeleton to be at least 20 um long
Img_Mask_VesselLumen_Skel(cat(1,rp_Mask_VesselLumen_Skel(FlagID).PixelIdxList)) = 0;

Img_Mask_VesselLumen = ...
    logical(imfill(Img_Mask_VL_SmoothBoundary,find(Img_Mask_VesselLumen_Skel))); % If skeleton closes in loop

if VesselLumenOnlyOption == 'n' 
    Img_Mask_Flow = ...
        logical(Img_Mask_Flow.* ... 
        imdilate(~Img_Tissue,strel('disk',round(1.0*CellDiameter_px_Input))));
    rp_Mask_Flow = regionprops(Img_Mask_Flow,'PixelIdxList');

    % Smooth boundary of thresholded blob area
    [Img_Mask_F_SmoothBoundary] = ...
        fSmoothImgBWBoundary_SG(Img_Mask_Flow,2,2*round(0.5*3.0*CellDiameter_px_Input)+1,XY_PxLength_Input); 
    rp_Mask_F_SmoothBoundary = regionprops(Img_Mask_F_SmoothBoundary,'PixelIdxList');
    rp_Mask_F_SmoothBoundary_PixelIdxList = cat(1,rp_Mask_F_SmoothBoundary.PixelIdxList);

    % Remove blobs from mask skipped by fSmoothImgBWBoundary_SG()
    FlagID = zeros(height(rp_Mask_Flow),1,'logical');
    for BlobCt = 1:height(rp_Mask_Flow)
        FlagID(BlobCt) = ~logical(numel(find(ismember(rp_Mask_F_SmoothBoundary_PixelIdxList,rp_Mask_Flow(BlobCt).PixelIdxList))));
    end
    Img_Mask_Flow(cat(1,rp_Mask_Flow(FlagID).PixelIdxList)) = 0;

    % Fill region within smoothed boundary; update Img_Mask_Flow
    Img_Mask_Flow = ...
        imfill(Img_Mask_F_SmoothBoundary,find(bwmorph(bwskel(Img_Mask_Flow,'MinBranchLength',round(4.0*CellDiameter_px_Input)),'spur',3)));
    
    % Fill residual holes in mask, if exist
    Img_Mask_Flow = imfill(Img_Mask_Flow,'holes');

    % Img_Mask_Flow must lie within Img_Mask_VesselLumen
    Img_Mask_Flow = logical(Img_Mask_Flow.*Img_Mask_VesselLumen);
else 
    Img_Mask_Flow = zeros(height(ImgStack_Input),width(ImgStack_Input),'logical');
end


%%

% Check
% % % figure; 
% % % Img_Temp = imfuse(Img_Mask_VesselLumen,Img_Mask_Flow,'Scaling','none','ColorChannels',[1 2 0]);
% % % imshow(Img_Temp,'InitialMagnification',100);












