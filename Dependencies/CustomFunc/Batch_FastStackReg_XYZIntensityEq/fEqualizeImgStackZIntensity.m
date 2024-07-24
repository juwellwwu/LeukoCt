function [ImgStack_Zeq_Output,Img_Mask_VesselLumen] = fEqualizeImgStackZIntensity(ImgStack_Input,XY_PxLength,BypassEqualizeZIntensityOption_Input)

if BypassEqualizeZIntensityOption_Input == 'n'

    fprintf('Measure and correct mean intensity difference between Z planes in stack...\n');
    
    [~,Img_Mask_VesselLumen,~] = fCreateRegImgZthresh(ImgStack_Input,XY_PxLength,'','y','n','n'); % default: 'y','n','n'

    % close all;
    % figure;
    % imshow(Img_Mask_VesselLumen);

    % == Set vessel areas to 0; only non-vessel areas of image compared
    ImgStack_NV = ImgStack_Input;
    ImgStack_NV(repmat(Img_Mask_VesselLumen,1,1,size(ImgStack_NV,3))) = 0;

    % == Remove close to saturated (or to 0) pixels
    NV_SatMask = ones(height(ImgStack_NV),width(ImgStack_NV),'logical');
    NV_SatMask(max(ImgStack_NV,[],3)>0.95) = 0;
    NV_SatMask(max(ImgStack_NV,[],3)<0.05) = 0;
    ImgStack_NV = ImgStack_NV.*repmat(NV_SatMask,1,1,size(ImgStack_NV,3));

    % == Sum intensity of each Z plane; normalize (set max = 1)
    NV_ZSum = sum(sum(ImgStack_NV,1),2);
    NV_ZSum = NV_ZSum./max(NV_ZSum,[],'all');

    % Plot Intensity change between Z planes
    % figure;
    % plot(permute(ImgStack_NV_ZPlaneSumNorm,[3 1 2]));

    % == Correct Z Sum intensity difference
    ImgStack_Zeq_Output = ImgStack_Input./repmat(NV_ZSum,size(ImgStack_Input,1),size(ImgStack_Input,2),1);
    ImgStack_Zeq_Output = ImgStack_Zeq_Output./max(ImgStack_Zeq_Output,[],'all');

else

    ImgStack_Zeq_Output = ImgStack_Input;
    Img_Mask_VesselLumen = zeros(height(ImgStack_Input),width(ImgStack_Input),'logical');
end


