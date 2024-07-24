clc; close all hidden; fclose('all');
% clearvars;


%% Save list of variables in Workspace after previous module, if does not exist

if ~exist('varList_MOD01') && ~exist('varList_MOD02') 
    varList_MOD01 = who;
    varList_MOD01 = varList_MOD01(~ismember(varList_MOD01,varList_All));
    varList_All = vertcat(varList_All,varList_MOD01);
end

if ~exist('DEMOOption')
    DEMOOption = 'n';
end


%% Create Output Folder

SaveMODFilePath = strcat(SaveFilePath,'STEP03_Velocity/');
if exist(SaveMODFilePath,'dir')==7
    rmdir(SaveMODFilePath,'s'); % important
end
mkdir(SaveMODFilePath);

if DEMOOption == 'y'
    SaveMODDEMOFilePath = strcat(SaveMODFilePath,'DEMOIllus/');
    mkdir(SaveMODDEMOFilePath);
    SkelCt_DEMO = 1;
    SegCt_DEMO = 3;
    BlkCt_DEMO = 1;
    gaussSigma_DEMO = 5.0;
    k_DEMO = 320;
end


%%  Blood Flow Velocity Measurement

for iCt = 1:3 % Iteration
    
    if iCt == 1
        clearvars OrientBandCt_gaussSigma_Cell;
        NumSkel_it = 1; 
        gaussSigma_Lo = 0.5;
        gaussSigma_Hi = (2*ceil(0.5*(0.5*min(oRadius.L.tMidPt_SkPx_Diameter,[],'all')/2.355))+1);
        gaussSigma_Step = (2.0:0.5:3.0)'; 
        [~,gaussSigma_Step_pdIdx] = pdist2(gaussSigma_Step,(gaussSigma_Hi-gaussSigma_Lo)/5,'euclidean','Smallest',1); 
        gaussSigma_Step = gaussSigma_Step(gaussSigma_Step_pdIdx);
        gaussSigma = (gaussSigma_Lo:gaussSigma_Step:gaussSigma_Hi)'; 
        WBC_MinVelocityEst_pxfr = 1500/XY_PxLength*Z_PxLength; % 1500 um/s Initial min WBC velocity guess
    elseif iCt == 2 
        OrientBandCt_gaussSigma_it1_Cell = OrientBandCt_gaussSigma_Cell;
        NumSkel_it = 1; 
        gaussSigma_Step = (1:0.5:2.0)'; 
        [~,gaussSigma_Step_pdIdx] = pdist2(gaussSigma_Step,(gaussSigma_Hi-gaussSigma_Lo)/8,'euclidean','Smallest',1); 
        gaussSigma_Step = gaussSigma_Step(gaussSigma_Step_pdIdx);
        gaussSigma = (gaussSigma_Lo:gaussSigma_Step:gaussSigma_Hi)';
        WBC_MinVelocityEst_pxfr = min(cell2mat(Velocity_Cell),[],'all','omitmissing'); 
    elseif iCt == 3
        OrientBandCt_gaussSigma_it2_Cell = OrientBandCt_gaussSigma_Cell; 
        NumSkel_it = NumSkel; 
        gaussSigma= (max((min(cell2mat(OrientBandCt_gaussSigma_it2_Cell),[],'all')-1),1):0.5:(max(cell2mat(OrientBandCt_gaussSigma_it2_Cell),[],'all')+1))'; 
        WBC_MinVelocityEst_pxfr = min(cell2mat(Velocity_Cell),[],'all','omitmissing'); 
    end


    %% Prepare "window" for suppressing replication artefacts during Fourier Transform

    Img_DR_FFT_Window = fCreateImgFFTWindow([TimeParam.tSeg_FrLength SkelBlockLength]); % [TimeParam.tSeg_FrLength SkelBlockLength] = size(Img_DR_Cell{1,1,1})

    if DEMOOption == 'y' 
        if (iCt==3) 
            imwrite(double(Img_DR_FFT_Window),strcat(SaveMODDEMOFilePath,'Img_DR_FFT_Window.tif'));
        end
    end


    %% Prepare frequency domain size mask 
    
    % = Low frequency cutoff
    if iCt ==1 
        CellDiameter_px_LoFreqCutoff = (12:2.0:18.0)./XY_PxLength; 
    elseif iCt ==2 
        CellDiameter_px_LoFreqCutoff = (max(1,floor(12.0/max(cell2mat(OrientBandCt_gaussSigma_it1_Cell),[],'all'))):2.0:18)./XY_PxLength; 
    elseif iCt ==3
        CellDiameter_px_LoFreqCutoff = (max(1,floor(12.0/max(cell2mat(OrientBandCt_gaussSigma_it2_Cell),[],'all'))):2.0:18)./XY_PxLength;
    end
    Img_DR_FFT_Mask_LoFreqCutoff = zeros(size(CellDiameter_px_LoFreqCutoff),'single');
    for cCt = 1:numel(CellDiameter_px_LoFreqCutoff) 
        BandWidth_fr = (min(CellDiameter_px_LoFreqCutoff(1,cCt)/WBC_MinVelocityEst_pxfr,0.8*min(oRadius.L.tMidPt_SkPx_Diameter,[],'all')/WBC_MinVelocityEst_pxfr));
        Img_DR_FFT_Mask_LoFreqCutoff(1,cCt) = round(TimeParam.tSeg_FrLength/BandWidth_fr);
        if Img_DR_FFT_Mask_LoFreqCutoff(1,cCt)>TimeParam.tSeg_FrLength/5 
            Img_DR_FFT_Mask_LoFreqCutoff(1,cCt) = TimeParam.tSeg_FrLength/5;
        end
    end

    [Img_DR_FFT_Mask_LoFreqCutoff, Img_DR_FFT_Mask_LoFreqCutoff_UniqueValIdx] = unique(Img_DR_FFT_Mask_LoFreqCutoff(1,cCt));
    CellDiameter_px_LoFreqCutoff = CellDiameter_px_LoFreqCutoff(Img_DR_FFT_Mask_LoFreqCutoff_UniqueValIdx);
    
    % = High frequency cutoff
    Img_DR_FFT_Mask_HiFreqCutoff = ... 
        round(SkelBlockLength/(1.00./XY_PxLength));

    Img_DR_FFT_FreqCutoff_MaskStack = fCreateFFTFreqCutoffMask(Img_DR_FFT_Mask_LoFreqCutoff,Img_DR_FFT_Mask_HiFreqCutoff,[TimeParam.tSeg_FrLength SkelBlockLength],2);

    if DEMOOption == 'y' 
        if (iCt==3)
            save(strcat(SaveMODDEMOFilePath,'Img_DR_FFT_FreqCutoff_MaskStack.mat'),'Img_DR_FFT_FreqCutoff_MaskStack');
            imwrite(double(Img_DR_FFT_FreqCutoff_MaskStack(:,:,1)),strcat(SaveMODDEMOFilePath,'Img_DR_FFT_FreqCutoff_MaskStack_Slice01.tif'));
        end
    end

    clearvars CellDiameter_px_LoFreqCutoff Img_DR_FFT_Mask_LoFreqCutoff_UniqueValIdx BandWidth_fr;


    %% Begin iteration by gaussSigma

    Img_DR_Postprocess_Cell_gsCell = cell(numel(gaussSigma),1);
    Img_DR_FFT_Filtered_Cell_gsCell = cell(numel(gaussSigma),1);
    Img_DR_FFT_Filtered_Clean_Cell_gsCell = cell(numel(gaussSigma),1);
    rp_FFT_Filtered_Clean_Cell_gsCell = cell(numel(gaussSigma),1);
    NumOrientBand_FFT_Filtered_Clean_Cell_gsCell = cell(numel(gaussSigma),1);

    for gsCt = 1:numel(gaussSigma) % # Gaussian sigma values

        %% Preprocess ImgStack for Creating Space-Time Diagram

        fprintf(strcat('\ngaussSigma =',32,num2str(gaussSigma(gsCt,1),'%0.2f'),'...\n'));

        ImgStack_Preprocess_Temp = ImgStack_Crop;
        ImgStack_Preprocess_Temp = imgaussfilt(ImgStack_Preprocess_Temp,gaussSigma(gsCt));
        ImgStack_Preprocess_Temp = ...
            (ImgStack_Preprocess_Temp-min(ImgStack_Preprocess_Temp,[],'all')).*(1-0)./(max(ImgStack_Preprocess_Temp,[],'all')-min(ImgStack_Preprocess_Temp,[],'all'));
        ImgStack_Preprocess_Temp = ImgStack_Preprocess_Temp.*oSCS.ImgStack_lThreshROI;

        % % % figure;
        % % % for k = 1:size(ImgStack_Crop,3)
        % % %     Img_Temp = ImgStack_Preprocess_Crop(:,:,k);
        % % %     % Img_Temp = imfuse(ImgStack_Crop(:,:,k),ImgStack_Preprocess_Crop(:,:,k),'ColorChannels',[1 2 0]);
        % % %     imshow(Img_Temp,'InitialMagnification',200);
        % % % end

        if DEMOOption == 'y'
            if (iCt==3) & abs(gaussSigma(gsCt)-gaussSigma_DEMO)<2.5
                imwrite(double(ImgStack_Preprocess_Temp(:,:,k_DEMO)),strcat(SaveMODDEMOFilePath,'ImgStack_Preprocess_Temp_iCt3_gaussSigma',num2str(gaussSigma(gsCt)),'_Slice',num2str(k_DEMO),'.tif'));
            end
        end

        if DEMOOption == 'y' 
            if (iCt==3) & (gsCt==1)
                imwrite(double(ImgStack_Crop(:,:,k_DEMO).*oSCS.ImgStack_lThreshROI(:,:,k_DEMO)),strcat(SaveMODDEMOFilePath,'ImgStack_Preprocess_Temp_iCt3_NoGaussFilt_Slice',num2str(k_DEMO),'.tif'));
            end
        end

        clearvars Img_Temp k;


        %% Create Space-Time Diagram

        Img_DR_Temp_Cell = ...
            fCreateImgDRCell(ImgStack_Preprocess_Temp,oSCS.Skel_xy_Cell,NumSkel_it,NumSeg,NumBlk,TimeParam);

        if DEMOOption == 'y' 
            if (iCt==3) & abs(gaussSigma(gsCt)-gaussSigma_DEMO)<2.5
                imwrite(double(Img_DR_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO}),strcat(SaveMODDEMOFilePath,'Img_DR_Temp_Cell_iCt3_gaussSigma',num2str(gaussSigma(gsCt)),'_Skel',num2str(SkelCt_DEMO,'%02d'),'_Seg',num2str(SegCt_DEMO,'%04d'),'_Blk',num2str(BlkCt_DEMO,'%02d'),'.tif'));
            end
        end

        if DEMOOption == 'y' 
            if (iCt==3) & (gsCt==1)
                Img_DR_NoGauss_Temp_Cell = ... % Generate Img_DR_Cell w/out Gaussian Filtering for Comparison
                    fCreateImgDRCell(ImgStack_Crop.*oSCS.ImgStack_lThreshROI,oSCS.Skel_xy_Cell,NumSkel_it,NumSeg,NumBlk,TimeParam);
                imwrite(double(Img_DR_NoGauss_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO}),strcat(SaveMODDEMOFilePath,'Img_DR_Temp_Cell_iCt3_NoGaussFilt_Skel',num2str(SkelCt_DEMO,'%02d'),'_Seg',num2str(SegCt_DEMO,'%04d'),'_Blk',num2str(BlkCt_DEMO,'%02d'),'.tif'));
                clearvars Img_DR_NoGauss_Temp_Cell;
            end
        end

        clearvars ImgStack_Preprocess_Temp;


        %% Postprocess Space-Time Diagram

        Img_DR_Postprocess_Temp_Cell =  Img_DR_Temp_Cell;

        % = Remove vertical background striations    
        Img_DR_Bkgd_Temp_Cell = cell(NumSkel_it,NumSeg,NumBlk);
        for SegCt = 1:NumSeg
            for SkelCt = 1:NumSkel_it
                for BlkCt = 1:NumBlk
                    Img_DR_Bkgd_Temp_Cell{SkelCt,SegCt,BlkCt} = mean(Img_DR_Postprocess_Temp_Cell{SkelCt,SegCt,BlkCt},1,'omitnan') + 1E-9; % 1E-9 to suppress any 0 val
                    Img_DR_Postprocess_Temp_Cell{SkelCt,SegCt,BlkCt} = Img_DR_Postprocess_Temp_Cell{SkelCt,SegCt,BlkCt}./Img_DR_Bkgd_Temp_Cell{SkelCt,SegCt,BlkCt};
                    Img_DR_Postprocess_Temp_Cell{SkelCt,SegCt,BlkCt} = Img_DR_Postprocess_Temp_Cell{SkelCt,SegCt,BlkCt}.*0.20; % set av of rows value = 0.20 (visibility)
                end
            end
        end

        % = Adjust intensity
        LoIntensity = min(cellfun(@(x) min(x,[],'all'),Img_DR_Postprocess_Temp_Cell,'UniformOutput',true),[],'all');
        HiIntensity = max(cellfun(@(x) max(x,[],'all'),Img_DR_Postprocess_Temp_Cell,'UniformOutput',true),[],'all');
        Img_DR_Postprocess_Temp_Cell = cellfun(@(x) (x-LoIntensity).*(1-0)./(HiIntensity-LoIntensity),Img_DR_Postprocess_Temp_Cell,'UniformOutput',false);

        if DEMOOption == 'y'
            if (iCt==3) & abs(gaussSigma(gsCt)-gaussSigma_DEMO)<2.5
                imwrite(double(Img_DR_Postprocess_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO}),strcat(SaveMODDEMOFilePath,'Img_DR_Postprocess_Temp_Cell_iCt3_gaussSigma',num2str(gaussSigma(gsCt)),'_Skel',num2str(SkelCt_DEMO,'%02d'),'_Seg',num2str(SegCt_DEMO,'%04d'),'_Blk',num2str(BlkCt_DEMO,'%02d'),'.tif'));
            end
        end

        clearvars Img_DR_Bkgd_Temp_Cell LoIntensity HiIntensity;
        clearvars SegCt SkelCt BlkCt;


        %% Perform Fourier Transform 

        fprintf('Performing Fourier Transform...\n');

        % Initiate cells with zero images
        Img_DR_FFT_Temp_Cell = repmat({zeros(size(Img_DR_Postprocess_Temp_Cell{1,1,1}),'single')},NumSkel_it,NumSeg,NumBlk);
        Img_DR_FFT_Filtered_Temp_Cell = repmat({zeros(size(Img_DR_Postprocess_Temp_Cell{1,1,1}),'single')},NumSkel_it,NumSeg,NumBlk);

        Img_DR_FFT_Vis_PreMask_Temp_Cell = repmat({zeros(size(Img_DR_Postprocess_Temp_Cell{1,1,1}),'single')},NumSkel_it,NumSeg,NumBlk);
        Img_DR_FFT_Vis_PostMask_Temp_Cell = repmat({zeros(size(Img_DR_Postprocess_Temp_Cell{1,1,1}),'single')},NumSkel_it,NumSeg,NumBlk);
        Img_DR_FFT_MaskStack_Temp_Cell = repmat({zeros(size(Img_DR_Postprocess_Temp_Cell{1,1,1}),'logical')},NumSkel_it,NumSeg,NumBlk); 
     
        for SegCt = 1:NumSeg
            for SkelCt = 1:NumSkel_it
                for BlkCt = 1:NumBlk
                        
                    if (iCt<3) | ((iCt==3) & any(ismember((OrientBandCt_gaussSigma_it1_Cell{1,SegCt,BlkCt}-gaussSigma_Step):0.5:(OrientBandCt_gaussSigma_it1_Cell{1,SegCt,BlkCt}+gaussSigma_Step),gaussSigma(gsCt))))
                    
                        % = Perform Fast Fourier Transform (FFT)               
                        Img_DR_FFT_Temp_Cell{SkelCt,SegCt,BlkCt} = fftshift(fft2(Img_DR_Postprocess_Temp_Cell{SkelCt,SegCt,BlkCt}.*Img_DR_FFT_Window));

                        % = Create 2 additional frequency domain masks
                        % Mask 2) Delete frequencies corresponding to near-vertical striations in xDyT image (ie, horizontal lines in FFT image)
                        % Mask 3) Clean up low intensity frequencies not relevant to velocity determination.

                        % Mask 2
                        Img_DR_FFT_Vis_ColSum_Mask_Temp = fCreateFFTColSumMask(Img_DR_FFT_Temp_Cell{SkelCt,SegCt,BlkCt});

                        % Mask 3
                        for cCt = 1:size(Img_DR_FFT_FreqCutoff_MaskStack,3) 

                            Img_DR_FFT_FreqCutoff_Mask_Temp = Img_DR_FFT_FreqCutoff_MaskStack(:,:,cCt);

                            Img_DR_FFT_Vis_FracMaxVal_Mask_Temp = fCreateFFTFracMaxValMask(Img_DR_FFT_Temp_Cell{SkelCt,SegCt,BlkCt},0.15);

                            % = Combine all masks into one stack
                            Img_DR_FFT_MaskStack_Temp = ...
                                cat(3,Img_DR_FFT_FreqCutoff_Mask_Temp,Img_DR_FFT_Vis_ColSum_Mask_Temp,Img_DR_FFT_Vis_FracMaxVal_Mask_Temp);

                            if sum(logical(min(Img_DR_FFT_MaskStack_Temp,[],3)),"all")>100
                                break 
                            end

                        end

                        Img_DR_FFT_Vis_PreMask_Temp_Cell{SkelCt,SegCt,BlkCt} = log10(abs(Img_DR_FFT_Temp_Cell{SkelCt,SegCt,BlkCt})); % For visualization
                        Img_DR_FFT_Vis_PreMask_Temp_Cell{SkelCt,SegCt,BlkCt}(Img_DR_FFT_Vis_PreMask_Temp_Cell{SkelCt,SegCt,BlkCt}<0) = 0;
                        Img_DR_FFT_MaskStack_Temp_Cell{SkelCt,SegCt,BlkCt} = Img_DR_FFT_MaskStack_Temp;
                        Img_DR_FFT_Vis_PostMask_Temp_Cell{SkelCt,SegCt,BlkCt} = ...  
                            Img_DR_FFT_Vis_PreMask_Temp_Cell{SkelCt,SegCt,BlkCt}.*min(Img_DR_FFT_MaskStack_Temp_Cell{SkelCt,SegCt,BlkCt},[],3);

                    end

                end
            end
        end

        if DEMOOption == 'y' 
            if (iCt==3) & abs(gaussSigma(gsCt)-gaussSigma_DEMO)<2.5
                imwrite(double(Img_DR_FFT_Vis_PreMask_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO}./max(Img_DR_FFT_Vis_PreMask_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO},[],'all')),strcat(SaveMODDEMOFilePath,'Img_DR_FFT_Vis_PreMask_Temp_Cell_iCt3_gaussSigma',num2str(gaussSigma(gsCt)),'_Skel',num2str(SkelCt_DEMO,'%02d'),'_Seg',num2str(SegCt_DEMO,'%04d'),'_Blk',num2str(BlkCt_DEMO,'%02d'),'.tif'));
                imwrite(double(Img_DR_FFT_Vis_PostMask_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO}./max(Img_DR_FFT_Vis_PostMask_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO},[],'all')),strcat(SaveMODDEMOFilePath,'Img_DR_FFT_Vis_PostMask_Temp_Cell_iCt3_gaussSigma',num2str(gaussSigma(gsCt)),'_Skel',num2str(SkelCt_DEMO,'%02d'),'_Seg',num2str(SegCt_DEMO,'%04d'),'_Blk',num2str(BlkCt_DEMO,'%02d'),'.tif'));
            end
        end

        clearvars cCt Img_DR_FFT_FreqCutoff_Mask_Temp Img_DR_FFT_Vis_ColSum_Mask_Temp Img_DR_FFT_Vis_FracMaxVal_Mask_Temp Img_DR_FFT_MaskStack_Temp;


        %% OPTIONAL: Save Img_DR_FFT_Vis_PreMask, Masks and Img_DR_FFT_Vis_PostMask for accuracy check

        if (iCt==3)
            if DEMOOption == 'n'
                fSaveFFTPrePostMask(Img_DR_FFT_Temp_Cell,Img_DR_FFT_MaskStack_Temp_Cell,'n','n',SaveMODFilePath);
            elseif DEMOOption == 'y'
                if (iCt==3) & (gaussSigma(gsCt)==gaussSigma_DEMO)
                    fSaveFFTPrePostMask(Img_DR_FFT_Temp_Cell,Img_DR_FFT_MaskStack_Temp_Cell,'n','y',SaveMODFilePath);
                end
            end
        end


        %% Return filtered Image to spatial domain by inverse Fourier Transform

        fprintf('Performing inverse Fourier Transform...\n');

        for SegCt = 1:NumSeg
            for SkelCt = 1:NumSkel_it
                for BlkCt = 1:NumBlk
                    
                    if (iCt<3) | ((iCt==3) & any(ismember((OrientBandCt_gaussSigma_it1_Cell{1,SegCt,BlkCt}-gaussSigma_Step):0.5:(OrientBandCt_gaussSigma_it1_Cell{1,SegCt,BlkCt}+gaussSigma_Step),gaussSigma(gsCt))))
                    
                        Img_DR_FFT_Temp_Cell{SkelCt,SegCt,BlkCt} = Img_DR_FFT_Temp_Cell{SkelCt,SegCt,BlkCt}.*logical(min(Img_DR_FFT_MaskStack_Temp_Cell{SkelCt,SegCt,BlkCt},[],3));
                        Img_DR_FFT_Filtered_Temp_Cell{SkelCt,SegCt,BlkCt} = ifft2(ifftshift(Img_DR_FFT_Temp_Cell{SkelCt,SegCt,BlkCt}), 'symmetric');

                    end

                end
            end
        end

        clearvars Img_FFT_Mask_LoFreqCutoff Img_FFT_Mask_HiFreqCutoff;
        clearvars Img_DR_FFT_Temp_Cell;


        %% Adjust brightness and contrast of Img_DR_FFT_Filtered_Cell
      
        Img_DR_FFT_Filtered_Temp_Cell = fAdjustImgFFTFilteredCellIntensity(Img_DR_FFT_Filtered_Temp_Cell,Img_DR_FFT_Vis_PreMask_Temp_Cell,Img_DR_FFT_Vis_PostMask_Temp_Cell);

         if DEMOOption == 'y' 
            if (iCt==3) & abs(gaussSigma(gsCt)-gaussSigma_DEMO)<2.5
                imwrite(double(Img_DR_FFT_Filtered_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO}),strcat(SaveMODDEMOFilePath,'Img_DR_FFT_Filtered_Temp_Cell_iCt3_gaussSigma',num2str(gaussSigma(gsCt)),'_Skel',num2str(SkelCt_DEMO,'%02d'),'_Seg',num2str(SegCt_DEMO,'%04d'),'_Blk',num2str(BlkCt_DEMO,'%02d'),'.tif'));
            end
         end


        %% Clean up Fourier Transform filtered Img_DR for band orientation and flow velocity measurement

        [Img_DR_FFT_Filtered_Clean_Temp_Cell,rp_FFT_Filtered_Clean_Temp_Cell] = ...
            fCleanImgFFTFilteredCell(Img_DR_FFT_Filtered_Temp_Cell,Img_DR_Postprocess_Temp_Cell,XY_PxLength);

         if DEMOOption == 'y' 
            if (iCt==3) & abs(gaussSigma(gsCt)-gaussSigma_DEMO)<2.5
                imwrite(double(Img_DR_FFT_Filtered_Clean_Temp_Cell{SkelCt_DEMO,SegCt_DEMO,BlkCt_DEMO}),strcat(SaveMODDEMOFilePath,'Img_DR_FFT_Filtered_Clean_Temp_Cell_iCt3_gaussSigma',num2str(gaussSigma(gsCt)),'_Skel',num2str(SkelCt_DEMO,'%02d'),'_Seg',num2str(SegCt_DEMO,'%04d'),'_Blk',num2str(BlkCt_DEMO,'%02d'),'.tif'));
            end
         end

        clearvars rp_FFT_Filtered_Clean_FlagBlobIdx_Cell;
        clearvars SkelCt SegCt MidptCt BlkCt rCt;


        %% Save different FFT filtered results based on different Gaussian smoothing of ImgStack_Crop

        Img_DR_Postprocess_Cell_gsCell{gsCt,1} = Img_DR_Postprocess_Temp_Cell;
        Img_DR_FFT_Filtered_Cell_gsCell{gsCt,1} = Img_DR_FFT_Filtered_Temp_Cell;
        Img_DR_FFT_Filtered_Clean_Cell_gsCell{gsCt,1} = Img_DR_FFT_Filtered_Clean_Temp_Cell;
        rp_FFT_Filtered_Clean_Cell_gsCell{gsCt,1} = rp_FFT_Filtered_Clean_Temp_Cell;
        NumOrientBand_FFT_Filtered_Clean_Cell_gsCell{gsCt,1} = cellfun(@height,rp_FFT_Filtered_Clean_Temp_Cell,'UniformOutput',true);


    end

    clearvars *_Temp*;


    %% Optimize gaussSigma for each SCS unit (SkelCt,SegCt,BlkCt) that maximizes the number of lines

    fprintf('Optimizing gaussSigma for each SCS unit that maximizes the number of lines...\n');

    NumOrientBand_FFT_Filtered_Clean_Cell_gsCell = cellfun(@(x) reshape(x,numel(x),1),NumOrientBand_FFT_Filtered_Clean_Cell_gsCell,'UniformOutput',false);

    [OrientBandCt,OrientBandCt_gaussSigma_maxIdx] = max(cell2mat(permute(NumOrientBand_FFT_Filtered_Clean_Cell_gsCell,[3,2,1])),[],3); % max by 3rd dim = gaussSigma

    NumOrientBand_FFT_Filtered_Clean_Cell_gsCell = cellfun(@(x) reshape(x,NumSkel_it,NumSeg,NumBlk),NumOrientBand_FFT_Filtered_Clean_Cell_gsCell,'UniformOutput',false);
    OrientBandCt_Cell = mat2cell(reshape(OrientBandCt,NumSkel_it,NumSeg,NumBlk),repelem(1,NumSkel_it),repelem(1,NumSeg),repelem(1,NumBlk));
    OrientBandCt_gaussSigma_maxIdx = reshape(OrientBandCt_gaussSigma_maxIdx,NumSkel_it,NumSeg,NumBlk);
    if isequal(size(OrientBandCt_gaussSigma_maxIdx),size(gaussSigma(OrientBandCt_gaussSigma_maxIdx)))
        OrientBandCt_gaussSigma_Cell = mat2cell(gaussSigma(OrientBandCt_gaussSigma_maxIdx),repelem(1,NumSkel_it),repelem(1,NumSeg),repelem(1,NumBlk));
    else % special case, when NumSkel_it==1 and NumBlk ==1
        OrientBandCt_gaussSigma_Cell = mat2cell(gaussSigma(OrientBandCt_gaussSigma_maxIdx)',repelem(1,NumSkel_it),repelem(1,NumSeg),repelem(1,NumBlk));
    end
    
    clearvars OrientBandCt;


    %% Build Img_DR_FFT_Filtered_Clean_Cell & rp_FFT_Filtered_Clean_Cell using optimal gaussSigma values

    if range(OrientBandCt_gaussSigma_maxIdx)<eps % same gaussSigma for all {SkelCt,SegCt,BlkCt}
        Img_DR_Postprocess_Cell = Img_DR_Postprocess_Cell_gsCell{OrientBandCt_gaussSigma_maxIdx(1,1,1),1}; 
        Img_DR_FFT_Filtered_Cell = Img_DR_FFT_Filtered_Cell_gsCell{OrientBandCt_gaussSigma_maxIdx(1,1,1),1}; 
        Img_DR_FFT_Filtered_Clean_Cell = Img_DR_FFT_Filtered_Clean_Cell_gsCell{OrientBandCt_gaussSigma_maxIdx(1,1,1),1};
        rp_FFT_Filtered_Clean_Cell = rp_FFT_Filtered_Clean_Cell_gsCell{OrientBandCt_gaussSigma_maxIdx(1,1,1),1};
    else
        for SegCt = 1:NumSeg
            for SkelCt = 1:NumSkel_it
                for BlkCt = 1:NumBlk

                    Img_DR_Postprocess_Cell{SkelCt,SegCt,BlkCt} = Img_DR_Postprocess_Cell_gsCell{OrientBandCt_gaussSigma_maxIdx(SkelCt,SegCt,BlkCt),1}{SkelCt,SegCt,BlkCt}; 
                    Img_DR_FFT_Filtered_Cell{SkelCt,SegCt,BlkCt} = Img_DR_FFT_Filtered_Cell_gsCell{OrientBandCt_gaussSigma_maxIdx(SkelCt,SegCt,BlkCt),1}{SkelCt,SegCt,BlkCt}; 

                    Img_DR_FFT_Filtered_Clean_Cell{SkelCt,SegCt,BlkCt} = Img_DR_FFT_Filtered_Clean_Cell_gsCell{OrientBandCt_gaussSigma_maxIdx(SkelCt,SegCt,BlkCt),1}{SkelCt,SegCt,BlkCt};
                    rp_FFT_Filtered_Clean_Cell{SkelCt,SegCt,BlkCt} = rp_FFT_Filtered_Clean_Cell_gsCell{OrientBandCt_gaussSigma_maxIdx(SkelCt,SegCt,BlkCt),1}{SkelCt,SegCt,BlkCt};

                end
            end
        end
    end

    clearvars gsCt SkelCt SegCt BlkCt;
    clearvars *_gsCell;


    %% Calculate Orient Band Orientation and Flow Velocity (px/fr)
    
    [oOrient, oVelocity] = ... 
        fCalcOrientVelocity(rp_FFT_Filtered_Clean_Cell,oSCS,TimeParam,[]); 
    Orient_Cell = oOrient.wMean_Cell; 
    Velocity_Cell = oVelocity.wMean_Cell; 


end

clearvars NumSkel_it gaussSigma_Lo gaussSigma_Hi gaussSigma_Step_pdIdx;
clearvars WBC_MinVelocityEst_pxfr gaussSigma gaussSigma_Step iCt OrientBandCt_gaussSigma_it1_Cell OrientBandCt_gaussSigma_it2_Cell;


%% OPTIONAL: Save Montage of gauss-Sigma optimized Orient Band slope lines  

SaveImgMontageFilePath = strcat(SaveMODFilePath,'Img_xDyT_OBline_BestgaussSigma/');
mkdir(SaveImgMontageFilePath);

LoHiIntensityLimit = stretchlim(reshape(cell2mat(Img_DR_FFT_Filtered_Cell),SkelBlockLength*TimeParam.tSeg_FrLength*numel(Img_DR_FFT_Filtered_Cell),1)); % based on all cells
Img_DR_FFT_Filtered_Cell_Save = cellfun(@(x) imadjust(x,LoHiIntensityLimit),Img_DR_FFT_Filtered_Cell,'UniformOutput',false);

fSaveImgDRXCellMontage(SaveImgMontageFilePath,Img_DR_Postprocess_Cell,Img_DR_FFT_Filtered_Cell_Save,Img_DR_FFT_Filtered_Clean_Cell);

clearvars Img_DR_FFT_Filtered_Cell;
clearvars LoHiIntensityLimit Img_DR_FFT_Filtered_Cell_Save SaveImgMontageFilePath;


%% Include gaussSigma Info in Orientation and Velocity output

oOrient.gaussSigma_Cell = OrientBandCt_gaussSigma_Cell;
oOrient.Img_DR_Postprocess_Cell = cellfun(@im2uint8,Img_DR_Postprocess_Cell,'UniformOutput',false); 
oOrient.Img_DR_OBline_Cell = Img_DR_FFT_Filtered_Clean_Cell; % for display
oOrient.rp_DR_OBline_Cell = rp_FFT_Filtered_Clean_Cell;


%% Clear variables likely no longer relevant after this module

clearvars Img_DR_FFT_Window;
clearvars Img_DR_FFT_Mask_LoFreqCutoff Img_DR_FFT_Mask_HiFreqCutoff Img_DR_FFT_FreqCutoff_MaskStack;
clearvars Orient_Cell OrientBandCt_Cell OrientBandCt_gaussSigma_Cell Img_DR_Postprocess_Cell Img_DR_FFT_Filtered_Clean_Cell rp_FFT_Filtered_Clean_Cell; 
clearvars Velocity_Cell; 
clearvars OrientBandCt_gaussSigma_maxIdx;
clearvars *_gsCell;
clearvars *DEMO*;


%% Save Workspace

% = OPTION 2: Save only new variables introduced
if ~exist('varList_MOD02') && ~exist('varList_MOD04') 
    varList_MOD02 = who;
    varList_MOD02 = varList_MOD02(~ismember(varList_MOD02,varList_All));
    varList_All = vertcat(varList_All,varList_MOD02);
    save(strcat(SaveMODFilePath,'Workspace_STEP03Var_',SampleIDString,'_',timestamp,'.mat'),varList_MOD02{:});
end


%% Save script in directory

ScriptName=mfilename;
PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveMODFilePath);
publish(strcat(ScriptName,'.m'),PublishOptions);

