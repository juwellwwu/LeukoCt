function [ImgStack_Reg_Output,RegShiftbyZ_Output,Img_Ref_Output] = ...
    fRegisterImgStack(ImgStack_Input,CropFlowMaskOption_Input,DoubleRegOption_Input,BypassRegOption_Input,XY_PxLength_Input)

ImgStack_Input = ImgStack_Input+(1E-6); 

stretchlim_Tol = [0.001 0.999];

Target_Height = 4*floor(0.25*height(ImgStack_Input));
Target_Width = 4*floor(0.25*width(ImgStack_Input));
cuboid = [floor(0.5*(width(ImgStack_Input)-Target_Width))+1 floor(0.5*(height(ImgStack_Input)-Target_Height))+1 1 Target_Width-1 Target_Height-1 size(ImgStack_Input,3)-1];
ImgStack_Input = imcrop3(ImgStack_Input,cuboid);

clearvars Target_Height Target_Width cuboid;


%% Registration
% Use File Exchange's "dftregistration" function for speed:
% https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

if BypassRegOption_Input == 'n'

    fprintf('Fast Registration of ImgStack_Input...\n');

    % == Elect Random Slice (k) to Register for 1st Pass
    k_rand_a = 0.25*size(ImgStack_Input,3);
    k_rand_b = 0.75*size(ImgStack_Input,3);
    k_rand = round((k_rand_b-k_rand_a).*rand(1,1)+k_rand_a);

    clearvars k_rand_a k_rand_b;

    % == Register
    fprintf('Registration starts.\n');

    if DoubleRegOption_Input == 'y'

        % = 1st pass (VERY ROUGH): do not register; calculate RegShiftbyZ so
        % know which k is optimal for registration

        fprintf('...REG: 1st pass...\n');

        ResizeScale = 0.25;
        ImgStack_Small = imresize(ImgStack_Input,ResizeScale,'Method','bilinear'); 
        
        Flatfield_sigma_Small = ResizeScale*(7.5)/XY_PxLength_Input/2.355; 
        Img_Ref_Small = ImgStack_Small(:,:,k_rand);
        Img_Ref_Small = imflatfield(Img_Ref_Small,Flatfield_sigma_Small);
        Img_Ref_Small = imadjust(Img_Ref_Small,stretchlim(Img_Ref_Small),stretchlim_Tol);

        RegShiftbyZ_Small = zeros(size(ImgStack_Input,3),4); 
        usfac = 2; 
        tic
        parfor k=1:size(ImgStack_Small,3) 
            ImgStack_Small(:,:,k) = imflatfield(ImgStack_Small(:,:,k),Flatfield_sigma_Small);
            ImgStack_Small(:,:,k) = imadjust(ImgStack_Small(:,:,k),stretchlim(ImgStack_Small(:,:,k)),stretchlim_Tol);
            [RegShiftbyZ_Small(k,:),~] = dftregistration_MOD(fft2(Img_Ref_Small),fft2(ImgStack_Small(:,:,k)),fft2(ImgStack_Small(:,:,k)),usfac);
        end
        toc

        clearvars ImgStack_Small Img_Ref_Small Flatfield_sigma_Small ResizeScale usfac k;

        % = 2nd Pass (ROUGH): 1st (ROUGH) Registration 
        % Use shift values determined w/ ImgStack_Small. This creates an almost registered stack quickly.

        fprintf('...REG: 2nd pass...\n');

        [~,k_MinShift] = min(sqrt((RegShiftbyZ_Small(:,3)-mean(RegShiftbyZ_Small(:,3))).^2 + (RegShiftbyZ_Small(:,4)-mean(RegShiftbyZ_Small(:,4))).^2));

        ResizeScale = 0.50;
        ImgStack_Small = imresize(ImgStack_Input,ResizeScale,'Method','bilinear'); 
        
        Flatfield_sigma_Small = ResizeScale*(7.5)/XY_PxLength_Input/2.355; 
        Img_Ref_Small = ImgStack_Small(:,:,k_MinShift);
        Img_Ref_Small = imflatfield(Img_Ref_Small,Flatfield_sigma_Small);
        Img_Ref_Small = imadjust(Img_Ref_Small,stretchlim(Img_Ref_Small),stretchlim_Tol);
        
        RegShiftbyZ_Small = zeros(size(ImgStack_Input,3),4);  
        usfac = 2*(1/ResizeScale); 
        tic
        parfor k=1:size(ImgStack_Small,3) 
            ImgStack_Small(:,:,k) = imflatfield(ImgStack_Small(:,:,k),Flatfield_sigma_Small);
            ImgStack_Small(:,:,k) = imadjust(ImgStack_Small(:,:,k),stretchlim(ImgStack_Small(:,:,k)),stretchlim_Tol);
            [RegShiftbyZ_Small(k,:),ImgStack_Reg_Output(:,:,k)] = ... 
                dftregistration_MOD(fft2(Img_Ref_Small),fft2(ImgStack_Small(:,:,k)),fft2(ImgStack_Input(:,:,k)),usfac);
            ImgStack_Reg_Output(:,:,k) = abs(ifft2(ImgStack_Reg_Output(:,:,k)));
        end
        toc
        
        % Crop area affected by registration
        [ImgStack_Reg_Output,~] = fCropImgStackReg(ImgStack_Reg_Output,RegShiftbyZ_Small,1,XY_PxLength_Input,CropFlowMaskOption_Input);

        clearvars ImgStack_Small Img_Ref_Small Flatfield_sigma_Small Flatfield_sigma_Small k_MinShift ResizeScale usfac k;

        % = 3rd pass, 2nd (FINE) Registration
        % Use reference Image k determined from 2nd pass

        fprintf('...REG: 3rd pass...\n');

        ImgStack_Reg_FF_Input = ImgStack_Reg_Output; 
        [gradThresh,numIter] = imdiffuseest(ImgStack_Reg_FF_Input(:,:,1),'ConductionMethod','quadratic');

        ResizeScale = 1.00;

        [~,k_MinShift] = min(sqrt((RegShiftbyZ_Small(:,3)-mean(RegShiftbyZ_Small(:,3))).^2 + (RegShiftbyZ_Small(:,4)-mean(RegShiftbyZ_Small(:,4))).^2));

        Flatfield_sigma = (7.5)/XY_PxLength_Input/2.355; 
        Img_Ref_Output = ImgStack_Reg_FF_Input(:,:,k_MinShift);
        Img_Ref_Output = imflatfield(Img_Ref_Output,Flatfield_sigma);
        Img_Ref_Output = imdiffusefilt(Img_Ref_Output,'ConductionMethod','quadratic', ...
            'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
        Img_Ref_Output = imadjust(Img_Ref_Output,stretchlim(Img_Ref_Output),stretchlim_Tol);

        RegShiftbyZ_Output = zeros(size(ImgStack_Reg_FF_Input,3),4); 
        usfac = 2; 
        tic
        parfor k=1:size(ImgStack_Reg_Output,3) 
            ImgStack_Reg_FF_Input(:,:,k) = ImgStack_Reg_Output(:,:,k);
            ImgStack_Reg_FF_Input(:,:,k) = imflatfield(ImgStack_Reg_FF_Input(:,:,k),Flatfield_sigma);
            ImgStack_Reg_FF_Input(:,:,k) = imdiffusefilt(ImgStack_Reg_FF_Input(:,:,k),'ConductionMethod','quadratic', ...
                'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
            ImgStack_Reg_FF_Input(:,:,k) = imadjust(ImgStack_Reg_FF_Input(:,:,k),stretchlim(ImgStack_Reg_FF_Input(:,:,k),stretchlim_Tol));
            [RegShiftbyZ_Output(k,:),ImgStack_Reg_Output(:,:,k)] = ... 
                dftregistration_MOD(fft2(Img_Ref_Output),fft2(ImgStack_Reg_FF_Input(:,:,k)),fft2(ImgStack_Reg_Output(:,:,k)),usfac);
            ImgStack_Reg_Output(:,:,k) = abs(ifft2(ImgStack_Reg_Output(:,:,k)));
        end
        toc

        clearvars ImgStack_Reg_FF_Input usfac k_MinShift ResizeScale Flatfield_sigma gradThresh numIter k;

        % Crop
        [ImgStack_Reg_Output,~] = fCropImgStackReg(ImgStack_Reg_Output,RegShiftbyZ_Output,1,XY_PxLength_Input,CropFlowMaskOption_Input);

    else

        % == 1st pass, 1st (FINE) Registration

        fprintf('...REG: 1st pass...\n');

        [gradThresh,numIter] = imdiffuseest(ImgStack_Input(:,:,1),'ConductionMethod','quadratic'); 
        
        Flatfield_sigma = (7.5)/XY_PxLength_Input/2.355; 
        Img_Ref_Output = ImgStack_Input(:,:,k_rand);
        Img_Ref_Output = imflatfield(Img_Ref_Output,Flatfield_sigma);
        Img_Ref_Output = imdiffusefilt(Img_Ref_Output,'ConductionMethod','quadratic', ...
                'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
        Img_Ref_Output = imadjust(Img_Ref_Output,stretchlim(Img_Ref_Output,stretchlim_Tol));

        RegShiftbyZ_Output = zeros(size(ImgStack_Input,3),4); 
        ImgStack_Reg_FF_Input = zeros(size(ImgStack_Input),'single'); 
        ImgStack_Reg_Output = zeros(size(ImgStack_Input),'single'); 
        usfac = 2; 
        tic
        parfor k=1:size(ImgStack_Reg_Output,3) % faster in parfor
            ImgStack_Reg_FF_Input(:,:,k) = ImgStack_Input(:,:,k);
            ImgStack_Reg_FF_Input(:,:,k) = imflatfield(ImgStack_Reg_FF_Input(:,:,k),Flatfield_sigma);
            ImgStack_Reg_FF_Input(:,:,k) = imdiffusefilt(ImgStack_Reg_FF_Input(:,:,k),'ConductionMethod','quadratic', ...
                'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
            ImgStack_Reg_FF_Input(:,:,k) = imadjust(ImgStack_Reg_FF_Input(:,:,k),stretchlim(ImgStack_Reg_FF_Input(:,:,k),stretchlim_Tol));
            [RegShiftbyZ_Output(k,:),~] = dftregistration_MOD(fft2(Img_Ref_Output),fft2(ImgStack_Reg_FF_Input(:,:,k)),fft2(ImgStack_Input(:,:,k)),usfac);
            ImgStack_Reg_Output(:,:,k) = abs(ifft2(ImgStack_Reg_Output(:,:,k)));
        end
        toc

        clearvars ImgStack_Reg_FF_Input Flatfield_sigma gradThresh numIter usfac k;
        
        % Crop
        [ImgStack_Reg_Output,~] = fCropImgStackReg(ImgStack_Reg_Output,RegShiftbyZ_Output,1,XY_PxLength_Input,CropFlowMaskOption_Input);

    end

else

    ImgStack_Reg_Output = ImgStack_Input;
    RegShiftbyZ_Output = nan(size(ImgStack_Input,3),4);
    Img_Ref_Output = zeros(size(ImgStack_Input(:,:,1)),'single');

end


%% Check (Optional)

% % % figure;
% % % for k = 1:size(ImgStack_Reg_Output,3)
% % %     imshow(ImgStack_Reg_Output(:,:,k));
% % %     % Img_Temp = imfuse(ImgStack_Input(:,:,k),ImgStack_Reg_Output(:,:,k),'ColorChannels',[1 2 0],'Scaling','none');
% % %     % imshow(Img_Temp,'InitialMagnification',200);
% % % end


%%

close all;

clearvars k_rand* k k_MinShift;
clearvars Img_Zstdev RectROI cuboid_RectROI;
clearvars ImgStack_RegAreaCrop;
clearvars Img_Zmin Img_ZmaxMask;
clearvars Img_Temp;