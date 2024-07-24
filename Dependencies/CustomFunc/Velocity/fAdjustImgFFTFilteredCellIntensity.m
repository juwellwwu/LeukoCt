function Img_FFT_Filtered_BCAdj_Cell_Output = fAdjustImgFFTFilteredCellIntensity(Img_FFT_Filtered_Cell_Input,Img_FFT_Vis_PreMask_Cell_Input,Img_FFT_Vis_PostMask_Cell_Input)

NumSkel_Input = size(Img_FFT_Filtered_Cell_Input,1);
NumSeg_Input = size(Img_FFT_Filtered_Cell_Input,2);
NumBlk_Input = size(Img_FFT_Filtered_Cell_Input,3);

% = Determine intensity scaling factor for each element
max_wrapper = @(x) max(x,[],'all','omitmissing');
Img_FFT_Vis_PreMask_MaxVal = cellfun(max_wrapper,Img_FFT_Vis_PreMask_Cell_Input);
Img_FFT_Vis_PreMask_MaxVal = Img_FFT_Vis_PreMask_MaxVal./max(Img_FFT_Vis_PreMask_MaxVal,[],'all');
Img_FFT_Vis_PostMask_MaxVal = cellfun(max_wrapper,Img_FFT_Vis_PostMask_Cell_Input);
Img_FFT_Vis_PostMask_MaxVal = Img_FFT_Vis_PostMask_MaxVal./max(Img_FFT_Vis_PostMask_MaxVal,[],'all');

Img_FFT_Vis_Scale = (1/max(Img_FFT_Vis_PostMask_MaxVal,[],"all",'omitmissing')).*(1./Img_FFT_Vis_PreMask_MaxVal);

% = Scale
Img_FFT_Filtered_BCAdj_Cell_Output = cell(NumSkel_Input,NumSeg_Input,NumBlk_Input);
for SegCt = 1:NumSeg_Input
    for SkelCt = 1:NumSkel_Input
        for BlkCt = 1:NumBlk_Input
            if Img_FFT_Vis_PreMask_MaxVal(SkelCt,SegCt,BlkCt)>eps 
                Img_FFT_Filtered_BCAdj_Cell_Output{SkelCt,SegCt,BlkCt} = Img_FFT_Filtered_Cell_Input{SkelCt,SegCt,BlkCt}.*Img_FFT_Vis_Scale(SkelCt,SegCt,BlkCt);
            else
                Img_FFT_Filtered_BCAdj_Cell_Output{SkelCt,SegCt,BlkCt} = zeros(size(Img_FFT_Filtered_Cell_Input{SkelCt,SegCt,BlkCt}),'single');
            end
        end
    end
end

% = Adjust max brightness such that max of Cell_Output = 1
MaxIntensity_Temp = max(cellfun(max_wrapper,Img_FFT_Filtered_BCAdj_Cell_Output),[],'all','omitmissing');
Img_FFT_Filtered_BCAdj_Cell_Output = ...
    cellfun(@(x) (x/MaxIntensity_Temp),Img_FFT_Filtered_BCAdj_Cell_Output,'UniformOutput',false);

% = Check outcome (Dim; need imadjust())
% figure;
% for SegCt = 1:NumSeg
%     for SkelCt = 1:NumSkel
%         for BlkCt = 1:NumBlk
%             imshow(Img_FFT_Filtered_BCAdj_Cell_Output{SkelCt,SegCt,BlkCt},'InitialMagnification',300);
%         end
%     end
% end

clearvars Img_DR_FFT_Vis_PreMask_MaxVal Img_DR_FFT_Vis_PostMask_MaxVal Img_DR_FFT_Vis_Scale;
clearvars max_wrapper MaxIntensity_Temp;