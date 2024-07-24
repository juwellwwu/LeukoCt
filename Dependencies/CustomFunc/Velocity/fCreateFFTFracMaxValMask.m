function [Img_FFT_FracMaxVal_Mask_Output] = fCreateFFTFracMaxValMask(Img_FFT,FracThreshold)

Img_FFT_PreMask = log10(abs(Img_FFT)); 
Img_FFT_PreMask(Img_FFT_PreMask<0) = 0; 

% = Create Mask
Img_FFT_FracMaxVal_Mask_Output = imgaussfilt(Img_FFT_PreMask,0.5); 
Img_FFT_FracMaxVal_Mask_Output = ...
    imbinarize(Img_FFT_FracMaxVal_Mask_Output./max(Img_FFT_FracMaxVal_Mask_Output,[],"all"),FracThreshold); 
Img_FFT_FracMaxVal_Mask_Output = bwmorph(Img_FFT_FracMaxVal_Mask_Output,'majority');

% = Gentle clean up
rp = regionprops(Img_FFT_FracMaxVal_Mask_Output,'Area','PixelIdxList');
FlagID = cat(1,rp.Area)<(0.05*max(cat(1,rp.Area),[],'all')) | cat(1,rp.Area)<4;
Img_FFT_FracMaxVal_Mask_Output(cat(1,rp(FlagID).PixelIdxList)) = 0;


