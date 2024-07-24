function [ImgStack_Out,XY_PxLength_Out,ZCropParam_Out,SampleIDString_Out] = ...
    fPrepTifVolume(TifVol_H5orTIF_FilenameString_In,XY_PxLength_In,AdjMagParam_In,ZCropParam_In,tSeg_FrLength_Input,expr_SampleID_In,CreateSampleIDStrOnlyOption)

if CreateSampleIDStrOnlyOption == 'n'

    %% Load Tif Volume

    fprintf('Loading Tif Volume...\n');

    [~,~,fExt] = fileparts(TifVol_H5orTIF_FilenameString_In);

    switch lower(fExt)
        case '.h5'
            % = HD5 Format
            ImgStack_Out = h5read(TifVol_H5orTIF_FilenameString_In,'/t0/channel0');
        case '.tif'
            % = TIF Single File Format
            ImgStack_Out = tiffreadVolume(TifVol_H5orTIF_FilenameString_In);
        otherwise
            error('Unexpected file extension: %s', fExt);
    end

    ImgStack_Out = im2single(ImgStack_Out);

    clearvars fExt;

    
    %% Update ZCropParam

    ZCropParam_Out = ZCropParam_In;
    if ZCropParam_Out.Option == 'n'
       ZCropParam_Out.ZStart = 1;
       ZCropParam_Out.ZEnd = tSeg_FrLength_Input*floor(size(ImgStack_Out,3)/tSeg_FrLength_Input);
       ZCropParam_Out.Option = 'y'; 
    else
       NumZSlices =  tSeg_FrLength_Input*floor((ZCropParam_Out.ZEnd-ZCropParam_Out.ZStart+1)/tSeg_FrLength_Input);
       ZCropParam_Out.ZEnd = ZCropParam_Out.ZStart+NumZSlices-1;
    end
    
    clearvars NumZSlices;

    
    %% Crop Tif Volume Z-slices based on ZCropParam input

    ImgStack_Out = ImgStack_Out(:,:,ZCropParam_Out.ZStart:ZCropParam_Out.ZEnd);


    %%  Rescale Tif Volume"s XY pixel size based on AdjMagParam input

    if AdjMagParam_In.Option == 'y'
        RescaleFactor = 1./(AdjMagParam_In.Target_XY_PxLength./XY_PxLength_In);
        ImgStack_Out = imresize(ImgStack_Out,RescaleFactor,'bilinear');
        XY_PxLength_Out = AdjMagParam_In.Target_XY_PxLength;
    else
        XY_PxLength_Out = XY_PxLength_In;
    end

    clearvars RescaleFactor;


    %% Adjust Brightness and Contrast

    ImgStack_Out = (ImgStack_Out-min(ImgStack_Out,[],'all')).*(1-0)./(max(ImgStack_Out,[],'all')-min(ImgStack_Out,[],'all'));

    [ImgStack_Out,~] = fEqualizeImgStackZIntensity(ImgStack_Out,XY_PxLength_In,'n');
    
    LoHigh = reshape(ImgStack_Out,numel(ImgStack_Out),1);
    LoHigh(LoHigh<(median(LoHigh,'all')-9*(1.4826*mad(LoHigh)))) = [];  
    LoHigh = stretchlim(LoHigh,[1E-5 (1-0)]);
    ImgStack_Out = (ImgStack_Out-LoHigh(1)).*(1-0)./(LoHigh(2)-LoHigh(1));

    
else

    ImgStack_Out = [];

end


%% Create SampleID string (for identification)

[StartIdx_SampleID,EndIdx_SampleID] = regexp(TifVol_H5orTIF_FilenameString_In,expr_SampleID_In);
if ~isempty(StartIdx_SampleID)
    SampleIDString_Out = TifVol_H5orTIF_FilenameString_In(StartIdx_SampleID:EndIdx_SampleID);
else
    SampleIDString_Out = '';
end

clearvars StartIdx_SampleID EndIdx_SampleID
