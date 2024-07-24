function [Img_FFT_ColSum_Mask_Output] = fCreateFFTColSumMask(Img_FFT)

Img_FFT_PreMask = log10(abs(Img_FFT));
Img_FFT_PreMask(Img_FFT_PreMask<0) = 0; 

Img_FFT_ColSum_Mask_Output = ...
    ~logical(repmat(isoutlier(sum(Img_FFT_PreMask,2),'median','ThresholdFactor',3.0),1,size(Img_FFT_PreMask,2)));