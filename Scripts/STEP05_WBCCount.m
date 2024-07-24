clc; close all hidden; fclose('all');
% clearvars;


%% Save list of variables in Workspace after previous module, if does not exist

if ~exist('varList_MOD04') && ~exist('varList_MOD05') % haven't run MOD4 yet
    varList_MOD04 = who;
    varList_MOD04 = varList_MOD04(~ismember(varList_MOD04,varList_All));
    varList_All = vertcat(varList_All,varList_MOD04);
end

if ~exist('DEMOOption')
    DEMOOption = 'n';
end


%% Elect Methods to process

% = GateXY_Method: "1" to run; '0" to Skip
% METHOD 4): FrIntensitySum based on best focused vessel lumen pixels in window
GateXY_Method = logical([0 0 0 1]); % Do not change

% = GateXY_SingleWinLength_Option
GateXY_SingleWinLength_Option = 0; % Do not change

% = GateXY_CellCt_Method
GateXY_CellCt_Method = 2; % Do not change


%% Supplement Information Not Entered During User Input (from Code Development)

if ~exist('WBCtoGateXYWinAreaRatio') | isempty(WBCtoGateXYWinAreaRatio_Init)
    WBCtoGateXYWinAreaRatio_Init = [0.40 0.35 0.30 0.25 0.20];
end

if ~exist('tFr_MeanTissueIntensity') | isempty(tFr_MeanTissueIntensity)
    [tFr_MeanTissueIntensity] = fCalctFrTissueIntensity(ImgStack_Crop,~oSCS.Img_lThreshROI,'y');
end

%  Load Manual Cell Ct Frame information
clearvars Fr_ManualCellCt;
ZCropParamMOD.Option = 'y'; % Fixed in code-

if ~exist('Fr_ManualCellCt') | isempty(Fr_ManualCellCt)
    [Fr_ManualCellCt] = fLoadManualCellCtFr(Fr_ManualCellCt_FoldernameString,SampleIDString,ZCropParamMOD);
end

clearvars *SampleID* -except SampleIDString;
clearvars fileID_Fr_ManualCellCt formatSpec_Fr_ManualCellCt size_Fr_ManualCellCt;
clearvars Fr_ManualCellCt_FilenameString;


%% Create Output Folder

SaveMODFilePath = strcat(SaveFilePath,'STEP05_WBCCount/');
if exist(SaveMODFilePath,'dir')==7
    rmdir(SaveMODFilePath,'s'); % important
end
mkdir(SaveMODFilePath);

if DEMOOption == 'y'
    SaveMODDEMOFilePath = strcat(SaveMODFilePath,'DEMOIllus/');
    mkdir(SaveMODDEMOFilePath);
    if ~exist('oROI.ImgStack_PreROICrop_Sz')
        oROI.ImgStack_PreROICrop_Sz = [221 208]; % Input if not present
    end
end


%% METHOD 4): FrIntensitySum based on best focused vessel lumen pixels in window

ImgStack_FrIntensitySum = ImgStack_Crop;

FrIntensitySum_mCell = cell(1,numel(GateXY_Method));
oFlag_RBCBkgd_Cell_mCell = cell(1,numel(GateXY_Method));

for mCt = 4:4 

    if mCt==4 
        NumBlkorWin = numel(WBCtoGateXYWinAreaRatio_Init);
        GateXY_ByBlk_Option = 0;
        FrIntensitySum_BlkorWin = fCreateFrIntensitySumStruct(ImgStack_FrIntensitySum,NumBlkorWin);
        Method_Str = 'PxWin';
        fprintf('METHOD 4): FrIntensitySum based on best focused vessel lumen pixels in window...\n');
    end

    if GateXY_Method(mCt)

        for itCt = 1:2

            if itCt == 1 % 1st iteration: Multiple WinArea, 1 WinLoc each tested
                WBCtoGateXYWinAreaRatio = WBCtoGateXYWinAreaRatio_Init;
                GateXY_SingleWinLength_Option = 0;
                GateXY_SingleWinLengthUserChoice_Option = 0; 
                GateXY_SingleWinLengthUserChoice_Idx = 0; 
                CtPkGMM_RunOption = 2; 
            elseif itCt == 2 % 2nd iteration: 1 WinArea, Multiple WinLoc tested
                WBCtoGateXYWinAreaRatio = WBCtoGateXYWinAreaRatio_Init;
                GateXY_SingleWinLength_Option = 1;
                GateXY_SingleWinLengthUserChoice_Option = 1; 
                GateXY_SingleWinLengthUserChoice_Idx = find(FrIntensitySum_BlkorWin.CellCt.BrightRank==1); 
                CtPkGMM_RunOption = 3;  
            end

            % = Locate best window
            [WinSkelPxIdx_Cell,~] = fCalcSkelPxIdxGateXYWinArea(ImgStack_FrIntensitySum,oSCS,WBCtoGateXYWinAreaRatio,XY_PxLength);

            [GateXY_Param,Fig_SkelPxStartIdx_BestWin] = ...
                fLocateOptimalGateXYWindow(ImgStack_FrIntensitySum,oFocusQ.Img_FocusDiameter_Idx,...
                WinSkelPxIdx_Cell,WBCtoGateXYWinAreaRatio,...
                oSCS,... 
                oRadius,...
                XY_PxLength,...
                GateXY_ByBlk_Option,... 
                GateXY_SingleWinLength_Option,... 
                GateXY_SingleWinLengthUserChoice_Option, GateXY_SingleWinLengthUserChoice_Idx);

            % % % if ~isempty(Fig_SkelPxStartIdx_BestWin)
            % % %     imwrite(Fig_SkelPxStartIdx_BestWin,strcat(SaveMODFilePath,'Fig_SkelPxStartIdx_BestWin_',timestamp,'.tif'));
            % % % end

            % = Update empty structure for data loading
            NumBlkorWin = width(GateXY_Param.WinArea_px); % Update
            FrIntensitySum_BlkorWin = fCreateFrIntensitySumStruct(ImgStack_FrIntensitySum,NumBlkorWin);

            % = Load FrIntensitySum data
            for BlkorWinCt = 1:NumBlkorWin

                if isnan(GateXY_Param.WinSkelPxStartIdx(1,BlkorWinCt)) % No optimal window
                    continue;
                end

                % = Create Mask
                if mCt==4 
                    [FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt},~] = fConvertSkelxytoMask(ImgStack_FrIntensitySum,GateXY_Param.WinSkelPxLinIdx_Cell{1,BlkorWinCt},'n');
                    FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt} = logical(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}.*oSCS.Img_fThreshROI);
                end

                if ~isempty(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt})
                    FrIntensitySum_BlkorWin.Img_Mask_Overlay_Cell{1,BlkorWinCt} = imfuse(ImgStack_FrIntensitySum(:,:,1),0.50*single(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}),'ColorChannels',[1 2 0],'Scaling','none');
                end

                % = Calculate FrIntensitySum
                FrIntensitySum_BlkorWin.IPRaw(:,BlkorWinCt) = permute(sum(sum(ImgStack_FrIntensitySum.*single(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}),1),2),[3 1 2]);
                FrIntensitySum_BlkorWin.IPRaw(:,BlkorWinCt) = single(FrIntensitySum_BlkorWin.IPRaw(:,BlkorWinCt)./sum(single(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}),'all'));

                % = Maximum pixel value in Dim3 (depth) of ImgStack_FrIntensitySum
                if size(ImgStack_FrIntensitySum,3) > 1/Z_PxLength*5
                    FrIntensityMax_BlkorWin_Temp = ...
                        movmax(ImgStack_FrIntensitySum.*single(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}),1/Z_PxLength*5,3,'Endpoints','discard'); 
                else 
                    FrIntensityMax_BlkorWin_Temp = mean(maxk(ImgStack_FrIntensitySum.*single(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}),10,3),3);
                end
                FrIntensityMax_BlkorWin_Temp = ... 
                    reshape(FrIntensityMax_BlkorWin_Temp,height(FrIntensityMax_BlkorWin_Temp)*width(FrIntensityMax_BlkorWin_Temp),1,size(FrIntensityMax_BlkorWin_Temp,3));
                FrIntensityMax_BlkorWin_Temp = sort(FrIntensityMax_BlkorWin_Temp,1,'descend'); 
                FrIntensityMax_BlkorWin_Temp = ... 
                    mean(FrIntensityMax_BlkorWin_Temp(1:ceil(0.10*sum(single(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}),'all')),:,:),1);
                FrIntensityMax_BlkorWin_Temp = cat(3, ... 
                    repmat(FrIntensityMax_BlkorWin_Temp(:,:,1),1,1,floor(0.5*abs(size(ImgStack_FrIntensitySum,3)-size(FrIntensityMax_BlkorWin_Temp,3)))),...
                    FrIntensityMax_BlkorWin_Temp,...
                    repmat(FrIntensityMax_BlkorWin_Temp(:,:,end),1,1,ceil(0.5*abs(size(ImgStack_FrIntensitySum,3)-size(FrIntensityMax_BlkorWin_Temp,3)))));
                
                FrIntensitySum_BlkorWin.IPZmax(:,BlkorWinCt) = permute(FrIntensityMax_BlkorWin_Temp,[3 2 1]);

                clearvars FrIntensityMax_BlkorWin_Temp;

                % = How many pixels "long" is mask
                FrIntensitySum_BlkorWin.MaskSkelLength_px(1,BlkorWinCt) = GateXY_Param.WinLength_px(1,BlkorWinCt); 
                FrIntensitySum_BlkorWin.MaskSkelLength_um(1,BlkorWinCt) = GateXY_Param.WinLength_px(1,BlkorWinCt).*XY_PxLength;

                % = Standard WBC (12.0um) area / Mask Area (approximate value)
                FrIntensitySum_BlkorWin.MaskWBCAreaFrac(1,BlkorWinCt) = GateXY_Param.WBCtoGateXYWinAreaRatio(1,BlkorWinCt);

                % = Vessel Diameter
                FrIntensitySum_BlkorWin.VesselDiameter_px(1,BlkorWinCt) = GateXY_Param.WinVesselDiameter_px(1,BlkorWinCt);
                FrIntensitySum_BlkorWin.VesselDiameter_um(1,BlkorWinCt) = FrIntensitySum_BlkorWin.VesselDiameter_px(1,BlkorWinCt).*XY_PxLength;

                % = FocusDiameterIdxMean
                FrIntensitySum_BlkorWin.FocusDiameterIdxMean(1,BlkorWinCt) = ...
                    sum(oFocusQ.Img_FocusDiameter_Idx.*single(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}),'all','omitmissing')./sum(single(FrIntensitySum_BlkorWin.Img_Mask_Cell{1,BlkorWinCt}),'all');

            end

            % = Perform Bkgd+Noise Removal + CellCounting

            if GateXY_CellCt_Method == 1 % [Removed]
            elseif GateXY_CellCt_Method == 2 % 

                % i. Background and Noise Removal
                multScMAD_CellCt = 3.5; 
                [FrIntensitySum_BlkorWin.IPBkgd,FrIntensitySum_BlkorWin.IPBkgdRmv,FrIntensitySum_BlkorWin.IPMAD,~] = ...
                    fAnalyzeFrIntensityProfile(FrIntensitySum_BlkorWin.IPRaw, ...
                    multScMAD_CellCt,...
                    GateXY_Param,TimeParam,XY_PxLength,Z_PxLength,oVelocity,Fr_ManualCellCt);

                % ii. Cell Counting
                [CellCtGMM_Struc,Fig_CellCtGMM_Cell,Fig_GMMfit_Cell] = ...
                    fCountFrIntensityProfilePeaksGMM(FrIntensitySum_BlkorWin.IPBkgd,FrIntensitySum_BlkorWin.IPBkgdRmv,FrIntensitySum_BlkorWin.IPMAD,...
                    FrIntensitySum_BlkorWin.IPRaw,... 
                    multScMAD_CellCt,... 
                    GateXY_Param,TimeParam,Z_PxLength, ...
                    oVelocity,...
                    Fr_ManualCellCt,...
                    CtPkGMM_RunOption);
                [FrIntensitySum_BlkorWin.CellCt] = CellCtGMM_Struc;
                CellCt_Method_Str = 'GMM';

                % iii. Flag locations of bright RBC bkgd [Removed]
                if itCt == 2

                    multScMAD_CellCt = 3.5; 
                    [oFlag_RBCBkgd_Cell,Fig_RBCBkgd_Flag_Cell] = ... % Dummy Flag
                        fCreateDUMMYRBCBkgdFlagCell(FrIntensitySum_BlkorWin.IPBkgdRmv,FrIntensitySum_BlkorWin.CellCt.IPMADCC,FrIntensitySum_BlkorWin.IPBkgd,...
                        multScMAD_CellCt,...
                        FrIntensitySum_BlkorWin.IPZmax,...
                        GateXY_Param,TimeParam,oSCS,'n');
                    oFlag_RBCBkgd_Cell_mCell{1,mCt} = oFlag_RBCBkgd_Cell;

                end

            end

            % = Save montage of mask
            if (mCt==4) & (itCt==2)
                imwrite(cell2mat(FrIntensitySum_BlkorWin.Img_Mask_Overlay_Cell),strcat(SaveMODFilePath,Method_Str,'_ImgMaskOverlay_WBCtoGateXYWinAreaRatio_',strjoin(cellfun(@(x) num2str(x),mat2cell(GateXY_Param.WBCtoGateXYWinAreaRatio,1,repelem(1,numel(GateXY_Param.WBCtoGateXYWinAreaRatio))),'UniformOutput',false),'_'),'_',timestamp,'.tif'));
            end
            
            % = Save mask for Demo Illustration
            if itCt==2
                if DEMOOption == 'y'
                    for BlkorWinCt = 1:NumBlkorWin
                        Img_Illus_PreROICrop_ImgMask = zeros(oROI.ImgStack_PreROICrop_Sz(1),oROI.ImgStack_PreROICrop_Sz(2),'logical'); % Size fits PreROICrop Image
                        Img_Illus_PreROICrop_ImgMask(oROI.ROICrop_cuboid(2):oROI.ROICrop_cuboid(2)+height(ImgStack_Crop)-1,...
                            oROI.ROICrop_cuboid(1):oROI.ROICrop_cuboid(1)+width(ImgStack_Crop)-1) = ...
                            FrIntensitySum_BlkorWin.Img_Mask_Overlay_Cell{1,BlkorWinCt}(:,:,2);
                        imwrite(uint8(255.*Img_Illus_PreROICrop_ImgMask),strcat(SaveMODDEMOFilePath,Method_Str,'_ImgMaskIllus_itCt2_BlkorWinCt',num2str(BlkorWinCt),'_',timestamp,'.tif'));
                        imwrite(uint8(255.*bwperim(Img_Illus_PreROICrop_ImgMask)),strcat(SaveMODDEMOFilePath,Method_Str,'_ImgMaskPerimIllus_itCt2_BlkorWinCt',num2str(BlkorWinCt),'_',timestamp,'.tif'));
                    end
                end
            end
            clearvars BlkorWinCt Img_Illus_PreROICrop_ImgMask;

            % = Keep Fr_ManualCellCt in structure, if not empty
            if ~isempty(Fr_ManualCellCt)
                FrIntensitySum_BlkorWin.Fr_ManualCellCt = Fr_ManualCellCt;
            end

        end

        % = Copy GateXY_Param into FrIntensitySum_BlkorWin.CellCt
        FrIntensitySum_BlkorWin.CellCt = catstruct(FrIntensitySum_BlkorWin.CellCt, GateXY_Param);

    end

    % = Specify FrIntensitySum_Blk as METHOD data
    FrIntensitySum_mCell{1,mCt} = FrIntensitySum_BlkorWin;

end

clearvars tMidPt_Velocity_EstWBCPkWidth_Cell;
clearvars Fig_SkelPxStartIdx_BestWin;
clearvars NumBlkorWin GateXY_ByBlk_Option FrIntensitySum_BlkorWin Method_Str mCt;
clearvars CellCtGMM_Struc Img_Fig_CellCtGMM_Cell Fig_CellCtGMM_Cell Img_Fig_GMMfit_Cell Fig_GMMfit_Cell;
clearvars Img_Fig_RBCBkgd_Flag_Cell Fig_RBCBkgd_Flag_Cell oFlag_RBCBkgd_Cell;
clearvars WBCtoGateXYWinAreaRatio GateXY_SingleWinLength_Option GateXY_SingleWinLengthUserChoice_Option GateXY_SingleWinLengthUserChoice_Idx CtPkGMM_RunOption itCt;
clearvars BlkorWinCt;


%% Save as Output structure

oCellCt = struct;

oCellCt.Ct_mCell = FrIntensitySum_mCell;
oCellCt.Method =  GateXY_Method; % 1 if run for 4 methods: SkelBlk, PxBlk, SkelWin, PxWin
oCellCt.Method_Str = {'SkelBlk','PxBlk','SkelWin','PxWin'};
oCellCt.CtMethod = GateXY_CellCt_Method; % "1": SG (sgolay filter); "2": GMM (Gaussian mixture model)

if exist('Fr_ManualCellCt')
    oCellCt.tFr_Manual = single(Fr_ManualCellCt); % list of tFr at which WBC is seen during manual count; size (# Manually-counted-WBC x 1)
    oCellCt.tFr_Manual_R_Ct = numel(oCellCt.tFr_Manual); % # WBC, w/ repeated tFr (such as doublets) counted as different cells; size (1 x 1 x 1)
    oCellCt.tFr_Manual_NR_Ct = numel(unique(oCellCt.tFr_Manual)); % # WBC, w/ repeated tFr (such as doublets) counted as one cell only; size (1 x 1 x 1)
else
    oCellCt.tFr_Manual = [];
    oCellCt.tFr_Manual_R_Ct = -1; % For reporting
    oCellCt.tFr_Manual_NR_Ct = -1; % For reporting
end

clearvars FrIntensitySum_mCell Fr_ManualCellCt;


%% Cell Ct Flag [Removed]

oFlag_CellCt_mCell = cell(1,numel(oCellCt.Method));

for mCt = 4:4

    if oCellCt.Method(1,mCt)

        oFlag_CellCt = struct; % Dummy

        if matches(oCellCt.Ct_mCell{1,mCt}.CellCt.Method(1,find(oCellCt.Ct_mCell{1,mCt}.CellCt.BrightRank==1)),"GMM") % if Best CellCt Window Ct Method is GMM
            oFlag_CellCt.tMidPt_Metric_Cell = ... % Same metric value for all cell elements
                repmat({single(-9999)},NumSkel_VsWidth,NumMidPt,NumBlk); % cell size (NumSkel_VsWidth x NumMidPt x NumBlk)
            oFlag_CellCt.tMidPt_1BlkAllSk_Metric_Cell = ... % Same metric value for all cell elements
                repmat({single(-9999)},1,NumMidPt,NumBlk); % cell size (1 x NumMidPt x NumBlk)
            oFlag_CellCt.tMidPtwMean_1BlkAllSk_Metric_Cell = ... % Same metric value for all cell elements
                repmat({single(-9999)},1,1,NumBlk); % cell size (1 x 1 x NumBlk)
        end

        oFlag_CellCt.tMidPt_1BlkAllSk_Flag_Cell = repmat({logical(0)},1,NumMidPt,NumBlk);
        oFlag_CellCt.tMidPt_1BlkAllSk_Flag_Ct = repmat({logical(0)},1,1,NumBlk);
        oFlag_CellCt.tMidPt_1BlkAllSk_Flag_Pct = repmat({single(0)},1,1,NumBlk);

        oFlag_CellCt.tFr_1BlkAllSk_Flag_Cell = repmat({repmat(logical(0),TimeParam.tMidPt_FrEnd(end),1)},1,1,NumBlk); % cell size (1 x 1 x NumBlk); each cell size (Numk x 1)
        oFlag_CellCt.tFr_1BlkAllSk_Flag_Ct = repmat({logical(0)},1,1,NumBlk);
        oFlag_CellCt.tFr_1BlkAllSk_Flag_Pct = repmat({single(0)},1,1,NumBlk);

        oFlag_CellCt.tFrMidPt_1BlkAllSk_Flag_Cell = oFlag_CellCt.tFr_1BlkAllSk_Flag_Cell;
        oFlag_CellCt.tFrMidPt_1BlkAllSk_Flag_Ct = oFlag_CellCt.tFr_1BlkAllSk_Flag_Ct;
        oFlag_CellCt.tFrMidPt_1BlkAllSk_Flag_Pct = oFlag_CellCt.tFr_1BlkAllSk_Flag_Pct;

        oFlag_CellCt.tMidPt_AllBlkAllSk_Flag_Ct = logical(0);
        oFlag_CellCt.tMidPt_AllBlkAllSk_Flag_Pct = single(0);

        oFlag_CellCt.tFr_AllBlkAllSk_Flag_Ct = logical(0);
        oFlag_CellCt.tFr_AllBlkAllSk_Flag_Pct = single(0);

        oFlag_CellCt.tFrMidPt_AllBlkAllSk_Flag_Ct = logical(0);
        oFlag_CellCt.tFrMidPt_AllBlkAllSk_Flag_Pct = single(0);

        oFlag_CellCt_mCell{1,mCt} = oFlag_CellCt;

    end

end

clearvars oFlag_CellCt Flag_val mCt BrightBlkorWinIdx;


%% Update tFr_Flag_Plot with Best Cell Ct Window's vessel block

tFr_Flag_Plot_mCell = cell(1,numel(GateXY_Method));

for mCt = 1:numel(GateXY_Method)

    if GateXY_Method(mCt)

        % Vessel block in which Best CellCt Window is Located
        WinBlk = oCellCt.Ct_mCell{1,mCt}.CellCt.WinSkelPxBlk(find(oCellCt.Ct_mCell{1,mCt}.CellCt.BrightRank==1));

        % Prepare struct for GateXY_Method
        tFr_Flag_Plot_mCell{1,mCt} = struct;

        % slanCM Colormap from File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap
        tFr_Flag_Plot_mCell{1,mCt}.cmap = colormap(slanCM('Pastel1',9)); % Match # Colors to create to # Colors in Swatch
        close all; % slanCM opens new figure;

        tFr_Flag_Plot_mCell{1,mCt}.Flag01.Data = oFlag_Radius.L.tFr_1BlkAllSk_Flag_Cell{1,1,WinBlk}; % Dummy
        tFr_Flag_Plot_mCell{1,mCt}.Flag01.Color = tFr_Flag_Plot_mCell{1,mCt}.cmap(1,:);

        tFr_Flag_Plot_mCell{1,mCt}.Flag02.Data = oFlag_VctyEr.OB.tFr_1BlkAllSk_Flag_Cell{1,1,WinBlk};
        tFr_Flag_Plot_mCell{1,mCt}.Flag02.Color = tFr_Flag_Plot_mCell{1,mCt}.cmap(2,:);

        tFr_Flag_Plot_mCell{1,mCt}.Flag03.Data = oFlag_VctyEr.SK.tFr_1BlkAllSk_Flag_Cell{1,1,WinBlk};
        tFr_Flag_Plot_mCell{1,mCt}.Flag03.Color = tFr_Flag_Plot_mCell{1,mCt}.cmap(3,:);

        tFr_Flag_Plot_mCell{1,mCt}.Flag04.Data = oFlag_FlowStb.tFr_1BlkAllSk_Flag_Cell{1,1,WinBlk}; % Dummy
        tFr_Flag_Plot_mCell{1,mCt}.Flag04.Color = tFr_Flag_Plot_mCell{1,mCt}.cmap(4,:);

        tFr_Flag_Plot_mCell{1,mCt}.Flag05.Data = ... % Dummy
            cell2mat(oFlag_RBCBkgd_Cell_mCell{1,mCt}{find(oCellCt.Ct_mCell{1,mCt}.CellCt.BrightRank)==1}.tFr_1BlkAllSk_Flag_Cell);
        tFr_Flag_Plot_mCell{1,mCt}.Flag05.Color = tFr_Flag_Plot_mCell{1,mCt}.cmap(5,:);

        tFr_Flag_Plot_mCell{1,mCt}.Flag06.Data = oFlag_CellCt_mCell{1,mCt}.tFr_1BlkAllSk_Flag_Cell{1,1,WinBlk}; % Dummy
        tFr_Flag_Plot_mCell{1,mCt}.Flag06.Color = tFr_Flag_Plot_mCell{1,mCt}.cmap(6,:);

    end

end

clearvars mCt WinBlk tFr_Flag_Prelim_Plot tFr_Flag_Plot;


%% Plot METHOD 4): FrIntensitySum based on best focused vessel lumen pixels in window

for mCt = 4:4 

    FrIntensitySum_BlkorWin = oCellCt.Ct_mCell{1,mCt};
    tFr_Flag_Plot = tFr_Flag_Plot_mCell{1,mCt};

    if mCt==4 && GateXY_Method(mCt) 
        fprintf('Plotting: METHOD 4) ...\n');
    elseif ~GateXY_Method(mCt)
        continue;
    end

    % Fig_a = fPlotFrIntensityProfile(FrIntensitySum_BlkorWin.IPBkgdRmv,FrIntensitySum_BlkorWin.CellCtRaw,oCellCt.tFr_Manual,tFr_Flag_Plot,3,1.000);
    Fig_b = fPlotFrIntensityProfile(FrIntensitySum_BlkorWin.CellCt.IPCC,FrIntensitySum_BlkorWin.CellCt,oCellCt.tFr_Manual,tFr_Flag_Plot,3,1.000); % Intensity Profile used for CellCounting

    if mCt==4
        % imwrite(Fig_a,strcat(SaveMODFilePath,oCellCt.Method_Str{1,mCt},'_',CellCt_Method_Str,'_Plot_IPBkgdRmv_WBCtoGateXYWinAreaRatio_',strjoin(cellfun(@(x) num2str(x),mat2cell(oCellCt.Ct_mCell{1,mCt}.CellCt.WBCtoGateXYWinAreaRatio,1,repelem(1,numel(oCellCt.Ct_mCell{1,mCt}.CellCt.WBCtoGateXYWinAreaRatio))),'UniformOutput',false),'_'),'_',timestamp,'.tif'));
        imwrite(Fig_b,strcat(SaveMODFilePath,oCellCt.Method_Str{1,mCt},'_',CellCt_Method_Str,'_Plot_IPCellCt_WBCtoGateXYWinAreaRatio_',strjoin(cellfun(@(x) num2str(x),mat2cell(oCellCt.Ct_mCell{1,mCt}.CellCt.WBCtoGateXYWinAreaRatio,1,repelem(1,numel(oCellCt.Ct_mCell{1,mCt}.CellCt.WBCtoGateXYWinAreaRatio))),'UniformOutput',false),'_'),'_',timestamp,'.tif'));
    end

end

close all hidden;
clearvars Fig_a Fig_b;
clearvars NumBlkorWin GateXY_ByBlk_Option FrIntensitySum_BlkorWin tFr_Flag_Plot Method_Str CellCt_Method_Str mCt;


%% Summarize, save cell count, vessel diameter (um) for vessel block or window, and FocalDiameterIdxSum for reporting

fprintf('Summarizing, saving cell count, vessel diameter (um) for vessel block or window, and FocalDiameterIdxSum for reporting ...\n');

oCellCt_Report_mCell = cell(1,numel(GateXY_Method));
oFlag_RBCBkgd_mCell = cell(1,numel(GateXY_Method));

for mCt = 1:numel(GateXY_Method)

    if GateXY_Method(mCt)

        oCellCt_Save = oCellCt.Ct_mCell{1,mCt};

        BrightBlkorWinIdx = find(oCellCt_Save.CellCt.BrightRank==1); 
        BlkSCS_BrightBlkorWinIdx = oCellCt_Save.CellCt.WinSkelPxBlk(find(oCellCt_Save.CellCt.BrightRank==1)); 

        oCellCt_Report_mCell{1,mCt} = {...
            oCellCt.Method_Str{1,mCt},... 
            oCellCt_Save.CellCt.Method(BrightBlkorWinIdx),... 
            oCellCt_Save.CellCt.CellCt(BrightBlkorWinIdx),... 
            oCellCt_Save.CellCt.CellCt(BrightBlkorWinIdx)/oVolume.V.tMidPtSum_AllSkOBwMean_nL(:,:,BlkSCS_BrightBlkorWinIdx),...  
            mean(oCellCt_Save.CellCt.pkZScore{BrightBlkorWinIdx},1),... 
            std(oCellCt_Save.CellCt.pkZScore{BrightBlkorWinIdx},0,1),... 
            prctile(oCellCt_Save.CellCt.pkZScore{BrightBlkorWinIdx},25),... 
            median(oCellCt_Save.CellCt.pkZScore{BrightBlkorWinIdx},1),... 
            prctile(oCellCt_Save.CellCt.pkZScore{BrightBlkorWinIdx},75),... 
            oCellCt_Save.CellCt.BrightRankVal(BrightBlkorWinIdx),... 
            oCellCt_Save.VesselDiameter_um(BrightBlkorWinIdx),...
            oCellCt_Save.MaskSkelLength_um(BrightBlkorWinIdx),...
            oCellCt_Save.MaskWBCAreaFrac(BrightBlkorWinIdx)};

        writecell(oCellCt_Report_mCell{1,mCt},strcat(SaveMODFilePath,'oCellCt_Summary_',oCellCt.Method_Str{1,mCt},'_',timestamp,'.csv'));

        oFlag_RBCBkgd_mCell{1,mCt} = oFlag_RBCBkgd_Cell_mCell{1,mCt}{1,BrightBlkorWinIdx};

    end

end

% = Incorporate oCellCt_Report_mCell into oCellCt for reporting
oCellCt.Report_mCell = oCellCt_Report_mCell;

clearvars BrightBlkorWinIdx BlkSCS_BrightBlkorWinIdx oCellCt_Save Method_Str mCt;
clearvars oFlag_RBCBkgd_Cell oFlag_RBCBkgd_Cell_mCell;
clearvars oCellCt_Report_mCell;


%% Clear variables not relevant after this module

clearvars multScMAD_CellCt;
clearvars GateXY_Param FrIntensitySum_mCell; 
clearvars WinSkelPxIdx_Cell;
clearvars ImgStack_FrIntensitySum ImgStack_Sm_Option;
clearvars GateXY_SingleWinLength_Option GateXY_Method GateXY_CellCt_Method;
clearvars tFr_Flag_Prelim_Plot;
clearvars tFr_MeanTissueIntensity;
clearvars DEMOOption SaveMODDEMOFilePath;
clearvars Blk_FocusDiameterIdxSum_MaxIdx BestBlkMethod BestBlk; 


%% Save list of variables in Workspace after previous module, if does not exist

% = OPTION 2: Save only new variables introduced
if ~exist('varList_MOD05') && ~exist('varList_MOD06') 
    varList_MOD05 = who;
    varList_MOD05 = varList_MOD05(~ismember(varList_MOD05,varList_All));
    varList_All = vertcat(varList_All,varList_MOD05);
    save(strcat(SaveMODFilePath,'Workspace_STEP05Var_',SampleIDString,'_',timestamp,'.mat'),varList_MOD05{:}); 
end


%% Save script in directory

ScriptName=mfilename;
PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveMODFilePath);
publish(strcat(ScriptName,'.m'),PublishOptions);

