function [Img_FFT_MaskStack_Output] = fCreateFFTFreqCutoffMask(Img_DR_FFT_Mask_LoFreqCutoff_Input,Img_DR_FFT_Mask_HiFreqCutoff_Input,size_Img,MaskOption)

if MaskOption==2

    Img_FFT_MaskStack_Output = ones(size_Img(1),size_Img(2),numel(Img_DR_FFT_Mask_LoFreqCutoff_Input),'logical');

    for rCt = 1:numel(Img_DR_FFT_Mask_LoFreqCutoff_Input)

        Img_DR_FFT_Mask = Img_FFT_MaskStack_Output(:,:,rCt);

        yCutoff = ...
            [(-(Img_DR_FFT_Mask_LoFreqCutoff_Input(1,rCt))+0.5*size_Img(1)) ((Img_DR_FFT_Mask_LoFreqCutoff_Input(1,rCt))+0.5*size_Img(1))];
        if all(ismember(yCutoff,1:size_Img(1))) 
            Img_DR_FFT_Mask(yCutoff(1):yCutoff(2),:) = 0;
        end

        xCutoff = ...
            [(-(Img_DR_FFT_Mask_HiFreqCutoff_Input)+0.5*size_Img(2)) ((Img_DR_FFT_Mask_HiFreqCutoff_Input)+0.5*size_Img(2))];
        if all(ismember(xCutoff,1:size_Img(2))) 
            Img_DR_FFT_Mask(:,xCutoff(1):xCutoff(2)) = 0;
        end

        Img_FFT_MaskStack_Output(:,:,rCt) = Img_DR_FFT_Mask;

    end

end

