function [ImgStack_XYeq_Output,Img_Mask_VesselLumen] = fEqualizeImgStackXYIntensity(ImgStack_Input,BkgdSmoothGaussSigma,XY_PxLength,BypassEqualizeXYIntensityOption_Input)

if BypassEqualizeXYIntensityOption_Input == 'n'

    fprintf('Measure and correct XY intensity difference in each Z plane of stack...\n');

    % == Create crude tissue mask with Z-project standard deviation

    [~,Img_Mask_VesselLumen,~] = fCreateRegImgZthresh(ImgStack_Input,XY_PxLength,'','y','n','n');

    Img_Mask_VesselLumen(:,1:2) = 0; 
    Img_Mask_VesselLumen(:,end-1:end) = 0;
    Img_Mask_VesselLumen(1:2,:) = 0;
    Img_Mask_VesselLumen(end-1:end,:) = 0;
    
    % close all;
    % figure;
    % imshow(Img_Mask_VesselLumen);

    % == Calculate ZMean in tissue area & correct XY intensity difference

    % = Create Vessel Lumen Masked Gaussian Background by slice
    ImgStack_XYeq_Bkgd = imresize(ImgStack_Input,0.5,'Method','bilinear'); 
    ImgStack_XYeq_Bkgd(repmat(imresize(Img_Mask_VesselLumen,0.5),1,1,size(ImgStack_Input,3))) = NaN;
    tic
    parfor k = 1:size(ImgStack_Input,3) 
        ImgStack_XYeq_Bkgd(:,:,k) = fillmissing2(ImgStack_XYeq_Bkgd(:,:,k),'linear'); % New function R2023a
        ImgStack_XYeq_Bkgd(:,:,k) = imgaussfilt(ImgStack_XYeq_Bkgd(:,:,k),0.5*BkgdSmoothGaussSigma);
    end
    toc
    
    % = Correct XY intensity difference
    ImgStack_XYeq_Output = ImgStack_Input./imresize(ImgStack_XYeq_Bkgd,[height(ImgStack_Input) width(ImgStack_Input)],'Method','bilinear');
    ImgStack_XYeq_Output = ImgStack_XYeq_Output./max(ImgStack_XYeq_Output,[],'all');

else 

    ImgStack_XYeq_Output = ImgStack_Input;
    Img_Mask_VesselLumen = zeros(height(ImgStack_Input),width(ImgStack_Input),'logical');

end
