clearvars;
clc; close all hidden; fclose('all');


%% User-Specified Variables
% Change to machine specifics to re

% === Input ImgStack folder: ImgStack can be h5 or TIF format *****
ImgStack_H5orTIF_Folder = '/Users/j/Documents/MATLAB/LeukoCt/Data_In/RAWVids_In/';

% === OBM Video Properties *****
VidProp_FilenameString = ...
    '/Users/j/Documents/MATLAB/LeukoCt/Data_In/OBMVideoProperties_In.xlsx';

% === Output image directory (do NOT include / at end) *****
OutputFolder = ...
    '/Users/j/Documents/MATLAB/LeukoCt/Data_Out/Reg/Batch_FastStackReg';

% ===  fRegisterImgStack_CropFlowMaskOption
fRegisterImgStack_CropFlowMaskOption = 'n'; % default: 'n'

% === Rotate
Rotate = 'n'; % default: 'n'
deg_Rotate = 0.0;

% === Adjust Magnification
AdjustMagnification_PreReg = 'n'; % default: 'n'
AdjustMagnification_PostReg = 'y'; % default: 'y'

% === Option to bypass registration, and/or background correction
BypassRegOption = 'n'; % default: 'n'
BypassEqualizeZIntensityOption = 'n'; % default: 'n'
BypassEqualizeXYIntensityOption = 'n'; % default: 'n'

% === Choose which stack to save
% 1: Save ImgStack_Save @ Reg_XY_PxLength Only
% 2: Save ImgStack_Save @ PostReg_XY_PxLength Only
% 3: Save ImgStack_Save @ Reg_XY_PxLength & PostReg_XY_PxLength
SaveOption = 3; % default: '3'


%% Create Output Folder for Batch

timestamp = datestr(datetime('now'),'yymmddHHMM');
SaveFilePath = strcat(OutputFolder,'_',timestamp,'/');
mkdir(SaveFilePath);


%% Identify Videos (=ImgStacks) in Batch Input Folder, load names

fprintf('Identifying ImgStack in batch input folder, loading names ...\n');

dirc = cat(1,dir(fullfile(ImgStack_H5orTIF_Folder,'*.h5')),dir(fullfile(ImgStack_H5orTIF_Folder,'*.tif')));
NumImgStack = height(dirc);
ImgStack_FilenameString_fCell = cell(1,NumImgStack);
Output_FoldernameString_fCell = cell(1,NumImgStack);
Output_FilenameString_fCell = cell(1,NumImgStack);

expr_SampleID = '(?-i)\d+_P\d+_V\d+_\d+(?:_S\d+-\d+)?'; % Sample ID, Optional Slice, ex 20231005_P15_V01_01_S01-99

% = "BatchRegSummary_fCell" contains reports of which run in batch input
% folder contains errors and is skipped, and at which task the error occurred
% Col 1: Sample ID
% Col 2: Task i-iv: 1 if Successful Run for SampleID, 0 otherwise
% * Task i: Registration
% * Task ii: Equalize intensity across each Z Plane (XY)
% * Task iii: Equalize intensity of Z Planes
% * Task iv: AdjustMagnification_PostReg
BatchRegSummary_fCell = repmat({"",zeros(1,4,'logical')},NumImgStack,1);

for fCt = 1:NumImgStack
    ImgStack_FilenameString_fCell{1,fCt} = fullfile(dirc(fCt,1).folder,dirc(fCt,1).name);
    [~,Output_FilenameString_fCell{1,fCt},~] = fileparts(dirc(fCt,1).name);
    Output_FoldernameString_fCell{1,fCt} = strcat(SaveFilePath,Output_FilenameString_fCell{1,fCt},'_FastStackReg_',timestamp,'/');
    mkdir(Output_FoldernameString_fCell{1,fCt});

    BatchRegSummary_fCell{fCt,1} = ImgStack_FilenameString_fCell{1,fCt};
end

clearvars dirc fCt;


%% Start Batch Registration

ImgStack_Output_Cell = cell(1,NumImgStack);
RegShiftbyZ_Cell = cell(1,NumImgStack);
Img_Ref_Cell = cell(1,NumImgStack);

fExt_Cell = cell(1,NumImgStack);

PreReg_XY_PxLength_Cell = cell(1,NumImgStack);
Reg_XY_PxLength_Cell = cell(1,NumImgStack);
PostReg_XY_PxLength_Cell = cell(1,NumImgStack);
XY_PxLength_Cell = cell(1,NumImgStack);
RescaleFactor_Cell = cell(1,NumImgStack);

Img_Mask_VesselLumen_EqXY_Cell = cell(1,NumImgStack);
Img_Mask_VesselLumen_EqZ_Cell = cell(1,NumImgStack);

for fCt = 1:NumImgStack

    tic

    %% Load individual ImgStack_Org

    fprintf(strcat('\nLoading ImgStack_Org',32,num2str(fCt,'%04d'),'/',num2str(NumImgStack,'%04d'),'...\n'))

    [~,~,fExt_Cell{1,fCt}] = fileparts(ImgStack_FilenameString_fCell{1,fCt});

    switch lower(fExt_Cell{1,fCt})
        case '.h5'
            % = HD5 Format
            ImgStack_Output_Cell{1,fCt} = h5read(ImgStack_FilenameString_fCell{1,fCt},'/t0/channel0');            
        case '.tif'
            % = TIF Format
            ImgStack_Output_Cell{1,fCt} = tiffreadVolume(ImgStack_FilenameString_fCell{1,fCt});
        otherwise
            error('Unexpected file extension: %s', fExt_Cell{1,fCt});
    end

    ImgStack_Output_Cell{1,fCt} = im2single(ImgStack_Output_Cell{1,fCt});


    %% Load Video Parameters: XY_PxLength, Z_PxLength
    % Load from Excel File

    fprintf('Loading pixel size information of ImgStack...\n');

    [~, ~, ~, PreReg_XY_PxLength_Cell{1,fCt}, Reg_XY_PxLength_Cell{1,fCt}, PostReg_XY_PxLength_Cell{1,fCt}] = ...
        fExtractVidProperty(VidProp_FilenameString, ImgStack_FilenameString_fCell{1,fCt}, expr_SampleID);


    %% OPTIONAL: Rotation

    if Rotate == 'y'
        [ImgStack_Output_Cell{1,fCt}] = fRotateCropImgStack(ImgStack_Output_Cell{1,fCt}, deg_Rotate);
    end


    %%  Rescale image; imresize resize 1st 2 dimensions of matrix

    if AdjustMagnification_PreReg == 'y'

        fprintf('Adjusting ImgStack"s magnification pre-registration...\n');

        RescaleFactor_Cell{1,fCt} = 1./(Reg_XY_PxLength_Cell{1,fCt}./PreReg_XY_PxLength_Cell{1,fCt});
        ImgStack_Output_Cell{1,fCt} = imresize(ImgStack_Output_Cell{1,fCt},RescaleFactor_Cell{1,fCt},'Method','bilinear');
        XY_PxLength_Cell{1,fCt} = Reg_XY_PxLength_Cell{1,fCt};

    else

        XY_PxLength_Cell{1,fCt} = PreReg_XY_PxLength_Cell{1,fCt};

    end

    ImgStack_Output_Cell{1,fCt} = ...
        (ImgStack_Output_Cell{1,fCt}-min(ImgStack_Output_Cell{1,fCt},[],'all')).*(1-0)./(max(ImgStack_Output_Cell{1,fCt},[],'all')-min(ImgStack_Output_Cell{1,fCt},[],'all'));


    %% Registration
    % Use File Exchange's "dftregistration" function:
    % https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

    try
        [ImgStack_Output_Cell{1,fCt},RegShiftbyZ_Cell{1,fCt},Img_Ref_Cell{1,fCt}] = ...
            fRegisterImgStack(ImgStack_Output_Cell{1,fCt},fRegisterImgStack_CropFlowMaskOption,'y',BypassRegOption,Reg_XY_PxLength_Cell{1,fCt});
        BatchRegSummary_fCell{fCt,2}(1,1) = 1; % Mark Success for Task i
    catch
        continue;
    end

    % % % ImgStack_Reg = ImgStack_Output_Cell{1,fCt};
    % % % RegShiftbyZ = RegShiftbyZ_Cell{1,fCt};


    %% Flat field correction: Equalize intensity across each Z plane (ie, in XY) of stack

    ImgStack_Output_Cell{1,fCt} = fCleanImgStackRegCrop(ImgStack_Output_Cell{1,fCt});

    try
        [ImgStack_Output_Cell{1,fCt},Img_Mask_VesselLumen_EqXY_Cell{1,fCt}] = ...
            fEqualizeImgStackXYIntensity(ImgStack_Output_Cell{1,fCt},45/XY_PxLength_Cell{1,fCt}/2.355,XY_PxLength_Cell{1,fCt},BypassEqualizeXYIntensityOption);
        BatchRegSummary_fCell{fCt,2}(1,2) = 1; % Mark Success for Task ii
    catch
        continue;
    end

    % % % ImgStack_EqXY = ImgStack_Output_Cell{1,fCt};


    %% Equalize intensity of Z planes in stack

    try
        [ImgStack_Output_Cell{1,fCt},Img_Mask_VesselLumen_EqZ_Cell{1,fCt}] = ...
            fEqualizeImgStackZIntensity(ImgStack_Output_Cell{1,fCt},XY_PxLength_Cell{1,fCt},BypassEqualizeZIntensityOption);
        BatchRegSummary_fCell{fCt,2}(1,3) = 1; % Mark Success for Task iii
    catch
        continue;
    end

    % % % ImgStack_EqZ = ImgStack_Output_Cell{1,fCt};


    %% Save translation data

    if ~(all(isnan(RegShiftbyZ_Cell{1,fCt})))
        writematrix(RegShiftbyZ_Cell{1,fCt},strcat(Output_FoldernameString_fCell{1,fCt},Output_FilenameString_fCell{1,fCt},'_RegShiftbyZ_',timestamp,'.csv'));
    end


    %% Adjust Intensity of final Image Stack, Save

    if (SaveOption==1 || SaveOption==3)

        fprintf('Adjusting Intensity of Final ImgStack; Saving...\n');

        % Final adjust intensity
        ImgStack_Output_Cell{1,fCt} = fAdjustImgStackIntensity(ImgStack_Output_Cell{1,fCt},XY_PxLength_Cell{1,fCt});

        % Save
        for k=1:size(ImgStack_Output_Cell{1,fCt},3)
            if k==1 
                imwrite(im2uint8(ImgStack_Output_Cell{1,fCt}(:,:,k)),strcat(Output_FoldernameString_fCell{1,fCt},Output_FilenameString_fCell{1,fCt},'_',num2str(1/XY_PxLength_Cell{1,fCt},'%0.2f'),'pxum_FastRegXYZIntEq.tif'));
            else
                imwrite(im2uint8(ImgStack_Output_Cell{1,fCt}(:,:,k)),strcat(Output_FoldernameString_fCell{1,fCt},Output_FilenameString_fCell{1,fCt},'_',num2str(1/XY_PxLength_Cell{1,fCt},'%0.2f'),'pxum_FastRegXYZIntEq.tif'),"WriteMode","append");
            end
        end

    end


    %% Post Registration Magnification Adjustment, Adjust Intensity and Save

    if AdjustMagnification_PostReg == 'y'

        try

            fprintf('Adjusting ImgStack"s magnification post-registration; Cropping and Saving ...\n');

            % Post Registration Magnification Adjustment
            RescaleFactor_Cell{1,fCt} = 1./(PostReg_XY_PxLength_Cell{1,fCt}./XY_PxLength_Cell{1,fCt});
            ImgStack_Output_Cell{1,fCt} = imresize(ImgStack_Output_Cell{1,fCt},RescaleFactor_Cell{1,fCt},'Method','bilinear');

            % Option: Crop to include vessel only, without too much extra space around.
            try
                ImgStack_Output_Cell{1,fCt} = fCropVesselArea(ImgStack_Output_Cell{1,fCt},PostReg_XY_PxLength_Cell{1,fCt});
                [ImgStack_Output_Cell{1,fCt},~] = fEqualizeImgStackZIntensity(ImgStack_Output_Cell{1,fCt},PostReg_XY_PxLength_Cell{1,fCt},'n'); % 20231024: Equalize Z intensity again after cropping
            catch
            end

            % Save as tiff volume
            if (SaveOption==2 || SaveOption==3)
                for k=1:size(ImgStack_Output_Cell{1,fCt},3)
                    if k==1
                        imwrite(im2uint8(ImgStack_Output_Cell{1,fCt}(:,:,k)),strcat(Output_FoldernameString_fCell{1,fCt},Output_FilenameString_fCell{1,fCt},'_',num2str(1/PostReg_XY_PxLength_Cell{1,fCt},'%0.2f'),'pxum_FastRegXYZIntEq.tif'));
                    else
                        imwrite(im2uint8(ImgStack_Output_Cell{1,fCt}(:,:,k)),strcat(Output_FoldernameString_fCell{1,fCt},Output_FilenameString_fCell{1,fCt},'_',num2str(1/PostReg_XY_PxLength_Cell{1,fCt},'%0.2f'),'pxum_FastRegXYZIntEq.tif'),"WriteMode","append");
                    end
                end
            end

            % Save as hd5
            if (SaveOption==2 || SaveOption==3)
                 h5create(strcat(Output_FoldernameString_fCell{1,fCt},Output_FilenameString_fCell{1,fCt},'_',num2str(1/PostReg_XY_PxLength_Cell{1,fCt},'%0.2f'),'pxum_FastRegXYZIntEq.h5'),...
                    '/t0/channel0',size(ImgStack_Output_Cell{1,fCt}),...
                    'Datatype','uint16',...
                    'ChunkSize',[height(ImgStack_Output_Cell{1,fCt}) width(ImgStack_Output_Cell{1,fCt}) max(1,floor(2^20/2/height(ImgStack_Output_Cell{1,fCt})/width(ImgStack_Output_Cell{1,fCt})))],...
                    'Deflate',1,'Shuffle',1,'Fletcher32',1); 
                h5write(strcat(Output_FoldernameString_fCell{1,fCt},Output_FilenameString_fCell{1,fCt},'_',num2str(1/PostReg_XY_PxLength_Cell{1,fCt},'%0.2f'),'pxum_FastRegXYZIntEq.h5'),'/t0/channel0',im2uint16(ImgStack_Output_Cell{1,fCt}));
                h5disp(strcat(Output_FoldernameString_fCell{1,fCt},Output_FilenameString_fCell{1,fCt},'_',num2str(1/PostReg_XY_PxLength_Cell{1,fCt},'%0.2f'),'pxum_FastRegXYZIntEq.h5')) % displays the metadata
            end

            BatchRegSummary_fCell{fCt,2}(1,4) = 1; 

        catch

            continue;

        end

    end

    % Clear from memory
    ImgStack_Output_Cell{1,fCt} = [];
    RegShiftbyZ_Cell{1,fCt} = [];
    pause(1);

    toc

end

clearvars fExt_Cell RescaleFactor_Cell XY_PxLength_Cell expr_SampleID;

% Save BatchRegSummary_fCell, which records which files ran successfully
writecell(BatchRegSummary_fCell,strcat(SaveFilePath,'BatchRegSummary_',timestamp,'.csv'));


%% Save script in directory
% html format; do not evaulate code or save figures

ScriptName=mfilename;
PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
publish(strcat(ScriptName,'.m'),PublishOptions);


%% =========== In-script Functions

function ImgStack_Out = fAdjustImgStackIntensity(ImgStack_In,XY_PxLength_In) 

LoHigh = stretchlim(reshape(ImgStack_In,numel(ImgStack_In),1),[1E-5 (1-1E-5)]);
ImgStack_Out = (ImgStack_In-LoHigh(1)).*(1-0)./(LoHigh(2)-LoHigh(1));

end

function ImgStack_Out = fCropVesselArea(ImgStack_In,XY_PxLength_In) 

[~,Img_Mask_VesselLumen,~] = fCreateRegImgZthresh(ImgStack_In,XY_PxLength_In,'','y','n','n');

rp_Mask = regionprops(Img_Mask_VesselLumen,'boundingbox');
rp_Mask_Bbox = cat(1,rp_Mask.BoundingBox);
rp_Mask_Bbox = ...
    horzcat(rp_Mask_Bbox,rp_Mask_Bbox(:,1)+rp_Mask_Bbox(:,3),rp_Mask_Bbox(:,2)+rp_Mask_Bbox(:,4));

rp_Mask_Bbox_Allxmin_Allymin = min(rp_Mask_Bbox(:,1:2),[],1); 
rp_Mask_Bbox_Allxmin_Allymin = rp_Mask_Bbox_Allxmin_Allymin-30; 
rp_Mask_Bbox_Allxmin_Allymin(rp_Mask_Bbox_Allxmin_Allymin<1) = 1;

rp_Mask_Bbox_Allxmax_Allymax = max(rp_Mask_Bbox(:,5:6),[],1); 
rp_Mask_Bbox_Allxmax_Allymax = rp_Mask_Bbox_Allxmax_Allymax+30; 
rp_Mask_Bbox_Allxmax_Allymax(1,1) = ... 
    min(rp_Mask_Bbox_Allxmax_Allymax(1,1),width(ImgStack_In));
rp_Mask_Bbox_Allxmax_Allymax(1,2) = ... 
    min(rp_Mask_Bbox_Allxmax_Allymax(1,2),height(ImgStack_In));

cuboid = [rp_Mask_Bbox_Allxmin_Allymin(1,1) rp_Mask_Bbox_Allxmin_Allymin(1,2) 1 (rp_Mask_Bbox_Allxmax_Allymax(1,1)-rp_Mask_Bbox_Allxmin_Allymin(1,1)) (rp_Mask_Bbox_Allxmax_Allymax(1,2)-rp_Mask_Bbox_Allxmin_Allymin(1,2)) size(ImgStack_In,3)-1];
ImgStack_Out = imcrop3(ImgStack_In,cuboid);

end

