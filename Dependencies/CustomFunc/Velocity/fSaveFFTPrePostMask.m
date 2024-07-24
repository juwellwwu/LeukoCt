function [] = fSaveFFTPrePostMask(Img_FFT_Cell,MaskStack_Cell,DisplayOption,SaveTIFOption,SaveMODFilePath_Input)

if DisplayOption=='n' && SaveTIFOption=='n'
    return 
end

if SaveTIFOption=='y'
    SaveFFTMaskFilePath = strcat(SaveMODFilePath_Input,'Img_xDyT_FFTMask/');
    if exist(SaveFFTMaskFilePath,'dir')==7
        rmdir(SaveFFTMaskFilePath,'s');
    end
    mkdir(SaveFFTMaskFilePath);
end

if DisplayOption=='y'
    figure;
end

NumSkel_Input = size(Img_FFT_Cell,1);
NumSeg_Input = size(Img_FFT_Cell,2);
NumBlk_Input = size(Img_FFT_Cell,3);

for SegCt = 1:NumSeg_Input
    for SkelCt = 1:NumSkel_Input
        for BlkCt = 1:NumBlk_Input

            Img_FFT_PreMask = log10(abs(Img_FFT_Cell{SkelCt,SegCt,BlkCt}));
            Img_FFT_PreMask(Img_FFT_PreMask<0) = 0; 

            % = All Mask
            Img_FFT_AllMask_Row = MaskStack_Cell{SkelCt,SegCt,BlkCt}(:,:,1); 
            for mCt = 2:size(MaskStack_Cell{SkelCt,SegCt,BlkCt},3) 
                Img_FFT_AllMask_Row = horzcat(Img_FFT_AllMask_Row,MaskStack_Cell{SkelCt,SegCt,BlkCt}(:,:,mCt));
            end
            Img_FFT_AllMask = logical(min(MaskStack_Cell{SkelCt,SegCt,BlkCt},[],3)); 

            % = Create image to display or save
            Img_FFT_PreMask = Img_FFT_PreMask./max(Img_FFT_PreMask,[],'all'); % max intensity = 1
            Img_FFT_PostMask = Img_FFT_PreMask.*Img_FFT_AllMask;
            Img_Temp = horzcat(repmat(Img_FFT_PreMask.*255,1,1,3),...
                repmat(Img_FFT_AllMask_Row.*255,1,1,3),...
                repmat(Img_FFT_PostMask.*255,1,1,3),imfuse(Img_FFT_PreMask,Img_FFT_PostMask,'ColorChannels',[1 2 0]));
            if DisplayOption=='y'
                imshow(Img_Temp,'InitialMagnification',300);
            end
            if SaveTIFOption=='y'
                imwrite(Img_Temp,strcat(SaveFFTMaskFilePath,'Img_DR_FFT_Pre_PostMask_Skel',num2str(SkelCt,'%02d'),'_Seg',num2str(SegCt,'%04d'),'_Blk',num2str(BlkCt,'%02d'),'.tif'));
            end

        end
    end
end

close all;

