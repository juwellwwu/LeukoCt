function [tFr_MeanTissueIntensity_Out] = fCalctFrTissueIntensity(ImgStack_In,TissueMask_In,TissueMaskCheck_Option)

if TissueMaskCheck_Option == 'y'
    Thresh_TissueIntensity_Temp = reshape(mean(ImgStack_In,3).*~TissueMask_In,numel(TissueMask_In),1);
    Thresh_TissueIntensity_Temp(Thresh_TissueIntensity_Temp<eps) = [];
    Thresh_TissueIntensity_Temp = prctile(Thresh_TissueIntensity_Temp,67); 
    TissueMask = imbinarize(mean(ImgStack_In,3),Thresh_TissueIntensity_Temp) & TissueMask_In;
else
    TissueMask = TissueMask_In;
end

maxk_k = round(0.50*sum(TissueMask,"all"));

ImgStack_Cell_In = mat2cell(ImgStack_In,height(ImgStack_In),width(ImgStack_In),repelem(1,size(ImgStack_In,3)));

% = Determine mean tissue intensity by z-plane 
tFr_MeanTissueIntensity_Out = cellfun(@(x) mean(maxk(reshape(x.*TissueMask,numel(TissueMask),1),maxk_k),'all'),ImgStack_Cell_In,'UniformOutput',false);
tFr_MeanTissueIntensity_Out = permute(cell2mat(tFr_MeanTissueIntensity_Out),[3 2 1]); 