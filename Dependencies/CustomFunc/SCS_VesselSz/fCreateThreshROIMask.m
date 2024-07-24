function [Img_ThreshROI_Output,ImgStack_ThreshROI_Output] = fCreateThreshROIMask(Img_Thresh_Cell_Input,Img_ROI_Input,TimeParam_Input)

Fr_Start_Input = TimeParam_Input.tSeg_FrStart;
Fr_End_Input = TimeParam_Input.tSeg_FrEnd;

ImgStack_Thresh_Input = cell2mat(permute(Img_Thresh_Cell_Input,[1,3,2]));

% == Time Independent ThreshROI Masks
Img_ThreshROI_Output = max(ImgStack_Thresh_Input,[],3).*Img_ROI_Input; 

% == Time-Dependent ThreshROI Mask
ImgStack_ThreshROI_Input = ImgStack_Thresh_Input.*repmat(Img_ROI_Input,1,1,size(ImgStack_Thresh_Input,3));

Numk = (size(ImgStack_ThreshROI_Input,3)+1)*(0.5*(Fr_End_Input(1)-Fr_Start_Input(1)+1)); 

ImgStack_ThreshROI_Output = zeros(height(ImgStack_ThreshROI_Input),width(ImgStack_ThreshROI_Input),Numk,'single');
for SegCt = 1:size(ImgStack_ThreshROI_Input,3)
    ImgStack_ThreshROI_Output(:,:,Fr_Start_Input(1,SegCt):Fr_End_Input(1,SegCt)) = ...
        ImgStack_ThreshROI_Output(:,:,Fr_Start_Input(1,SegCt):Fr_End_Input(1,SegCt)) + ImgStack_ThreshROI_Input(:,:,SegCt);
end

ImgStack_ThreshROI_Output(:,:,(Fr_Start_Input(1):Fr_Start_Input(2)-1)) = ... 
    ImgStack_ThreshROI_Output(:,:,(Fr_Start_Input(1):Fr_Start_Input(2)-1)).*2;

ImgStack_ThreshROI_Output(:,:,((Fr_End_Input(end-1)+1):Fr_End_Input(end))) = ... 
    ImgStack_ThreshROI_Output(:,:,((Fr_End_Input(end-1)+1):Fr_End_Input(end))).*2;

ImgStack_ThreshROI_Output(isnan(ImgStack_ThreshROI_Output)) = 0;
ImgStack_ThreshROI_Output(ImgStack_ThreshROI_Output<(1.0+eps)) = 0; 
ImgStack_ThreshROI_Output = logical(ImgStack_ThreshROI_Output);