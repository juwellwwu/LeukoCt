function [CellCtGMM_Struc_Out,Fig_CellCtGMM_Cell_Out,Fig_GMMfit_Cell_Out] = ...
    fCountFrIntensityProfilePeaksGMM(LinePlot_Bkgd_AllBlkorWin_In,LinePlot_BkgdRmv_AllBlkorWin_In,LinePlot_MAD_AllBlkorWin_In,...
    LinePlot_AllBlkorWin_In,... 
    multScMAD_CellCt_In,... 
    GateXYParam_In,TimeParam_In,Z_PxLength_In,...
    Velocity_Param_In,...
    Fr_ManualCellCt_In,... 
    CtPkGMM_RunOption_In) 

NumBlk_In = width(LinePlot_BkgdRmv_AllBlkorWin_In);

tMidPtSkelAv_Velocity_Cell_In = Velocity_Param_In.tMidPt_AllSkOBwMean_Cell; 


%% NumtFr For Smoothing

NumFr_HR = round((60/Z_PxLength_In)/240);


%% Multiple of time-dependent MAD used as cell count cutoff

multScMAD_CC = multScMAD_CellCt_In;


%% Set up Conditions of Run / Skip GMM depending on RunOption

if CtPkGMM_RunOption_In ==1 
    NumRun = 1;
    SkipGMM = zeros(1,NumBlk_In,'logical');
elseif CtPkGMM_RunOption_In == 2 
    NumRun = 1;
    SkipGMM = ones(1,NumBlk_In,'logical');
elseif CtPkGMM_RunOption_In == 3
    NumRun = 2;
    SkipGMM = ones(NumRun,NumBlk_In,'logical'); 
end

for rCt = 1:NumRun

    CellCt_AllBlkorWin = -9999*ones(1,width(LinePlot_AllBlkorWin_In),'single');
    CellCt_ExpEr = -9999*ones(1,width(LinePlot_AllBlkorWin_In),'single');
  
    CC_pkLoc_AllBlkorWin_Cell = cell(1,NumBlk_In); 
    CC_pkZScore_AllBlkorWin_Cell = cell(1,NumBlk_In); 
    CC_pktMidPt_AllBlkorWin_Cell = cell(1,NumBlk_In); 
    CC_pkParam_AllBlkorWin_Cell = cell(1,NumBlk_In); 

    CC_tMidPt_Ct_AllBlkorWin_Cell = cell(1,NumBlk_In);
    CC_rpGateXY_AllBlkorWin_Cell = cell(1,NumBlk_In);

    CC_IPCC_AllBlkorWin = -9999*ones(size(LinePlot_AllBlkorWin_In),'single');
    CC_IPMADCC_AllBlkorWin = -9999*ones(size(LinePlot_AllBlkorWin_In),'single');
    CC_IPBkgdRmvPkAg_AllBlkorWin = -9999*ones(size(LinePlot_AllBlkorWin_In),'single');

    Fig_CellCtGMM_Cell_Out = cell(1,NumBlk_In);
    Fig_GMMfit_Cell_Out = cell(1,NumBlk_In);

    for BlkCt = 1:NumBlk_In

        %% Load BlkorWin specific data
        LinePlot_In = LinePlot_AllBlkorWin_In(:,BlkCt);
        IP7_Bkgd = LinePlot_Bkgd_AllBlkorWin_In(:,BlkCt);
        IP7_BkgdRmv = LinePlot_BkgdRmv_AllBlkorWin_In(:,BlkCt);
        IP7_MAD = LinePlot_MAD_AllBlkorWin_In(:,BlkCt);

        WinLength_px_In = GateXYParam_In.WinLength_px(:,BlkCt);

        %% Expected WBC peak width in tFr

        tMidPtSkelAv_WinPkWidthEst_Cell = cellfun(@(x) WinLength_px_In/x,tMidPtSkelAv_Velocity_Cell_In,'UniformOutput',false);


        %% IP7_CellCt: Intensity Profile for Cell Counting

        IP7_MAD_CC = IP7_MAD; 

        % == Create IP7_CellCt: two-step creation
        % = 1) Create IP7_BkgdRmv_pkAug
        multScMAD_pkAug = multScMAD_CC*0.9;
        IP7_BkgdRmv_pkAug = IP7_BkgdRmv; 
        [pks,locs,w,p,wx] = findpeaks_MOD(IP7_BkgdRmv_pkAug,'WidthReference','halfheight');
        pkParam_Temp = horzcat(pks,locs,w,p,wx);
        FlagID = (pkParam_Temp(:,1)-multScMAD_pkAug*1.4826*IP7_MAD_CC(pkParam_Temp(:,2)))<0;
        pkParam_Temp(FlagID,:) = [];

        multScMAD_tr = 2.33; 
        [pks_tr,locs_tr,w_tr,p_tr,wx_tr] = findpeaks_MOD(-IP7_BkgdRmv_pkAug,'WidthReference','halfheight');
        trParam_Temp = horzcat(pks_tr,locs_tr,w_tr,p_tr,wx_tr);
        FlagID_tr = (trParam_Temp(:,1)-multScMAD_tr*1.4826*IP7_MAD_CC(trParam_Temp(:,2)))<0;
        trParam_Temp(FlagID_tr,:) = [];
        trParam_Temp(:,1) = -trParam_Temp(:,1);

        FlagID_tr = zeros(size(IP7_BkgdRmv_pkAug),'logical');

        for pkCt = 1:height(pkParam_Temp)

            LeftLim = pkParam_Temp(pkCt,2);
            RightLim = min(pkParam_Temp(pkCt,2)+8.0/1.1774*(ceil(pkParam_Temp(pkCt,6)-pkParam_Temp(pkCt,2))),height(LinePlot_In));

            if pkCt<height(pkParam_Temp)
                tr_pk = trParam_Temp(ismember(trParam_Temp(:,2),LeftLim:RightLim) & trParam_Temp(:,2)<pkParam_Temp(pkCt+1,2),:);
            else
                tr_pk = trParam_Temp(ismember(trParam_Temp(:,2),LeftLim:RightLim),:);
            end

            if ~isempty(tr_pk)

                [minVal,minIdx] = min(tr_pk(:,1),[],'all');
                pkValAugment = (abs(minVal)-multScMAD_tr*1.4826*IP7_MAD_CC(tr_pk(minIdx,2)));
                IP7_BkgdRmv_pkAug(floor(pkParam_Temp(pkCt,5)):ceil(pkParam_Temp(pkCt,6))) = IP7_BkgdRmv_pkAug(floor(pkParam_Temp(pkCt,5)):ceil(pkParam_Temp(pkCt,6)))+pkValAugment;

                LeftLim = max(1,tr_pk(minIdx,2)-2.0/1.1774*(tr_pk(minIdx,2)-tr_pk(minIdx,5)));
                RightLim = min(tr_pk(minIdx,2)+2.0/1.1774*(tr_pk(minIdx,6)-tr_pk(minIdx,2)),height(LinePlot_In));
                FlagID_Temp = find(IP7_BkgdRmv_pkAug<(-multScMAD_tr*1.4826*IP7_MAD_CC));
                FlagID_Temp(FlagID_Temp<LeftLim) = [];
                FlagID_Temp(FlagID_Temp>RightLim) = [];
                FlagID_tr(FlagID_Temp) = 1;

            end

        end

        IP7_BkgdRmv_pkAug(find(FlagID_tr)) = -multScMAD_tr*1.4826*IP7_MAD_CC(find(FlagID_tr));

        % 2) Create IP7_CellCt
        IP7_CellCt = IP7_BkgdRmv_pkAug-multScMAD_CC*1.4826*IP7_MAD_CC;
        IP7_CellCt(IP7_CellCt<0) = 0;
        
        close all;
        clearvars ax* tlo; 
        clearvars sgolay_FrLength sgolay_Order;
        clearvars pkCt tr_pk pkParam_Temp trParam_Temp minVal minIdx pkValAugment pks locs w p wx pks_tr locs_tr w_tr p_tr wx_tr LeftLim RightLim multScMAD_tr;


        %% Prepare for Gaussian Mixture Modeling: i) Quick Locate WBC peak regions w/ findpeaks()

        IP7_CellCt_FP = IP7_CellCt; 

        % == Prepare Crop Gate: keep all + values in IP7_CellCt_FP
        CropGate_CellCt_FP_IP7 = logical(IP7_CellCt_FP);
        rp_CropGate_CellCt_FP_IP7 = regionprops(CropGate_CellCt_FP_IP7,'Area','PixelIdxList');

        % == FPParam_IP7 will hold findpeaks() information
        % Col 1: peak height (pks)
        % Col 2: peak location (locs)
        % Col 3: peak width (w)
        % Col 4: peak prominence (p)
        % Col 5: peak width extent, Left Limit
        % Col 6: peak width extent, Right Limit
        % Col 7: peak Modified Z-score
        % Col 8: (Wlll load next step) peak location (locs) idx in tq_full scale
        % Col 9: peak location (locs)'s tMidPt

        FPParam_IP7 = -9999*ones(height(rp_CropGate_CellCt_FP_IP7)*3,9,'single');
        
        CellCtFP_IP7 = 0; 

        for gCt = 1:height(rp_CropGate_CellCt_FP_IP7) 

            [pks,locs,w,p,wx] = findpeaks_MOD(IP7_CellCt_FP,'WidthReference','halfheight');

            FlagID = ~ismember(locs,rp_CropGate_CellCt_FP_IP7(gCt).PixelIdxList);

            FPParam_IP7_g = horzcat(pks,locs,w,p,wx); 
            FPParam_IP7_g(FlagID,:) = [];

            FPParam_IP7_g = ...
                horzcat(FPParam_IP7_g,multScMAD_CC+FPParam_IP7_g(:,1)./(1.4826*IP7_MAD_CC(FPParam_IP7_g(:,2))));

            FPParam_IP7((CellCtFP_IP7+1):(CellCtFP_IP7+height(FPParam_IP7_g)),1:7) = FPParam_IP7_g;

            CellCtFP_IP7 = CellCtFP_IP7+height(FPParam_IP7_g);

        end

        FPParam_IP7(FPParam_IP7(:,1)<0,:) = [];
        FPParam_IP7(:,8) = 0;

        % = Include peak loc's tMidpt (Col 9)
        tMidPt_HistEdge = horzcat(TimeParam_In.tMidPt_FrStart,TimeParam_In.tMidPt_FrEnd(end));
        [tMidPt_HistCt,~,tMidPt_HistBin] = histcounts(FPParam_IP7(:,2),tMidPt_HistEdge);
        FPParam_IP7(:,9) = tMidPt_HistBin;

        tMidPt_IP7_CellCt_Cell = mat2cell(tMidPt_HistCt,1,repelem(1,numel(tMidPt_HistCt)),1); 

        % = Mean Z Score in each tMidPt
        if ~isempty(FPParam_IP7)
            [grp_tMidPt_HistBin,grpID_tMidPt_HistBin] = findgroups(tMidPt_HistBin);
            grp_ZScoreMean = splitapply(@mean, FPParam_IP7(:,7),grp_tMidPt_HistBin); 
            tMidPt_IP7_MeanZScore_Cell = zeros(1,TimeParam_In.Num_tMidPt,1);
            tMidPt_IP7_MeanZScore_Cell(grpID_tMidPt_HistBin) = grp_ZScoreMean;
            tMidPt_IP7_MeanZScore_Cell = mat2cell(tMidPt_IP7_MeanZScore_Cell,1,repelem(1,numel(tMidPt_HistCt),1));
        else
            tMidPt_IP7_MeanZScore_Cell = zeros(1,TimeParam_In.Num_tMidPt,1);
            tMidPt_IP7_MeanZScore_Cell = mat2cell(tMidPt_IP7_MeanZScore_Cell,1,repelem(1,numel(tMidPt_HistCt),1));
        end

        % Eliminate peaks after findpeaks() if peak height at FWHM
        % (left and right side) within 5% value of maximum peak height
        FlagID = (abs(floor(FPParam_IP7(:,5))-FPParam_IP7(:,2))>0) & ... 
            ((abs(IP7_CellCt_FP(floor(FPParam_IP7(:,5)))-IP7_CellCt_FP(FPParam_IP7(:,2)))./FPParam_IP7(:,1))<0.05);
        FPParam_IP7(FlagID,:) = [];

        FlagID = (abs(ceil(FPParam_IP7(:,6))-FPParam_IP7(:,2))>0) & ... 
            ((abs(IP7_CellCt_FP(ceil(FPParam_IP7(:,6)))-IP7_CellCt_FP(FPParam_IP7(:,2)))./FPParam_IP7(:,1))<0.05);
        FPParam_IP7(FlagID,:) = [];

        clearvars IP7_CellCt_FP; % = IP7_CellCt; not used again
        clearvars gCt CellCtFP_IP7 FPParam_IP7_g;
        clearvars tMidPt_Hist* grp*;

        if ~SkipGMM(rCt,BlkCt)


            %% Prepare for Gaussian Mixture Modeling: ii) Decrease time step size
            % New time step = tq_step, in unit of tFr

            IP7_CellCt_GMM = IP7_CellCt; 

            % = Make smaller time-step for analysis
            tq_step = 0.1;
            t_full = single(linspace(1,height(LinePlot_In),height(LinePlot_In))); 
            tq_full = ...
                single(linspace(1-(tq_step*floor(0.5*(1/tq_step))),height(LinePlot_In)+(tq_step*floor(0.5*(1/tq_step))),(1/tq_step)*height(LinePlot_In))); 
            CCq_full = interp1(t_full,IP7_CellCt_GMM,tq_full,'linear'); 
            CCq_full(CCq_full<(1E-3*max(CCq_full,[],'all'))) = 0; 
            CCq_full = single(fillmissing(CCq_full,'nearest')); 
            MDq_full = interp1(t_full,IP7_MAD_CC,tq_full,'linear'); 
            MDq_full = single(fillmissing(MDq_full,'nearest')); 

            % = Prepare CropGate in new time axis
            CropGate_CCq_full = logical(CCq_full); % All + values
            rp_CCq_full = regionprops(logical(CCq_full),CCq_full,'Area','PixelIdxList','PixelValues'); 
            FlagID = cat(1,rp_CCq_full.Area)<floor(0.5*(1/tq_step)); 
            rp_CCq_full(FlagID) = [];

            % = Determine FPParam_IP7 peak loc's index in new time axis
            [~,FPParam_IP7_locs_q] = pdist2(tq_full',FPParam_IP7(:,2),'euclidean','Smallest',1);
            FPParam_IP7(:,8) = FPParam_IP7_locs_q';
            
            clearvars FlagID;


            %% Gaussian Mixture Modeling (GMM), by each cropped video segment

            % GMMParam_IP7 stores Gaussian component properties and has 8 columns:
            % 1) gCt: which gate each Gaussian component belongs to
            % 2) GMM_mu: GMM model Gaussian component means
            % 3) GMM_sigma: GMM model Gaussian component sigmas
            % 4) GMM_weight: GMM model Gaussian component weights
            % 5) GMM_pdfMax: GMM model Gaussian's pdf function max ("peak") value
            % 6) GMM_pdfArea: GMM model Gaussian's pdf function area under curve
            % 7) GateArea: Area under gate's IP7_CellCt curve
            % 8) GMM_mu Z-Score

            GMMParam_IP7 = -9999*ones(height(FPParam_IP7).*5,8,'single');
            CellCt_GMM = 0; 
            Fig_GMMfit_gCell = cell(height(rp_CCq_full),1); 

            for gCt = 1:height(rp_CCq_full) 

                % === Prepare cropped, smaller-time-step intensity profile for GMM
                tq_g = tq_full(min(rp_CCq_full(gCt).PixelIdxList,[],'all'):max(rp_CCq_full(gCt).PixelIdxList,[],'all'));
                CCq_g = CCq_full(min(rp_CCq_full(gCt).PixelIdxList,[],'all'):max(rp_CCq_full(gCt).PixelIdxList,[],'all'));
                MDq_g = MDq_full(min(rp_CCq_full(gCt).PixelIdxList,[],'all'):max(rp_CCq_full(gCt).PixelIdxList,[],'all'));

                % === Calculate area under gate
                Area_g = trapz(tq_g,CCq_g); 

                % === Prepare input data for GMM
                histdata = [];
                for itCt = 1:2 % 
                    if itCt==1
                        multMDq = 1.0; 
                    elseif itCt==2
                        multMDq = 0.0;
                        CCq_g = ksdensity(histdata,tq_g,'Kernel','normal').*Area_g; 
                    end
                    data_multiple = round(500/prctile(CCq_g(CCq_g>0),10));
                    for tqCt = 1:numel(tq_g)
                        if CCq_g(tqCt)>0
                            histdata = ...
                                vertcat(histdata,(multMDq*1.4826*MDq_g(tqCt))*randn(round(data_multiple*CCq_g(tqCt)),1,'single')+tq_g(tqCt));
                        end
                    end
                end

                clearvars data_multiple multMDq tqCt itCt;

                % === Determine # Gaussian components (peaks) to fit within crop-gated stretch of intensity profile
                GMM_options = statset('Display','final','MaxIter',1000);

                % = Range of # Gaussian components to test
                [pks,locs,w,p,wx] = findpeaks_MOD(CCq_g','WidthReference','halfheight');
                pkParam_Temp = horzcat(pks,locs,w,p,wx);

                if height(pkParam_Temp)>1
                    FlagID = zeros(height(pkParam_Temp),1,'logical');
                    for pkCt = 1:height(pkParam_Temp)
                        if pkCt<2
                            MinPeakProminenceLim = 0.05*min([...
                                pkParam_Temp(pkCt,1)+multScMAD_CellCt_In*1.4826*MDq_g(pkParam_Temp(pkCt,2)),...
                                pkParam_Temp(pkCt+1,1)+multScMAD_CellCt_In*1.4826*MDq_g(pkParam_Temp(pkCt+1,2))]);
                        elseif pkCt>1 & pkCt<height(pkParam_Temp)
                            MinPeakProminenceLim = 0.05*min([...
                                pkParam_Temp(pkCt-1,1)+multScMAD_CellCt_In*1.4826*MDq_g(pkParam_Temp(pkCt-1,2)),...
                                pkParam_Temp(pkCt,1)+multScMAD_CellCt_In*1.4826*MDq_g(pkParam_Temp(pkCt,2)),...
                                pkParam_Temp(pkCt+1,1)+multScMAD_CellCt_In*1.4826*MDq_g(pkParam_Temp(pkCt+1,2))]);
                        elseif pkCt>(height(pkParam_Temp)-1)
                            MinPeakProminenceLim = 0.05*min([...
                                pkParam_Temp(pkCt-1,1)+multScMAD_CellCt_In*1.4826*MDq_g(pkParam_Temp(pkCt-1,2)),...
                                pkParam_Temp(pkCt,1)+multScMAD_CellCt_In*1.4826*MDq_g(pkParam_Temp(pkCt,2))]);
                        end
                        [pks,locs,w,p,wx] = findpeaks_MOD(CCq_g','WidthReference','halfheight','MinPeakProminence',MinPeakProminenceLim);
                        if ~isempty(locs)
                            FlagID = FlagID | ~ismember(pkParam_Temp(:,2),locs);
                        end
                    end
                    pkParam_Temp(FlagID,:) = [];
                end

                NumGMMComp_Est = height(pkParam_Temp);

                % = Prepare initial guess of Gaussian component for each
                % observation in histdata
                [~,fitgmdist_Start] = pdist2(tq_g(pkParam_Temp(:,2))',histdata,'euclidean','Smallest',1);
                for pkCt = 1:height(pkParam_Temp) 
                    LeftLim = tq_g(max(1,floor(pkParam_Temp(pkCt,5))));
                    RightLim = tq_g(min(ceil(pkParam_Temp(pkCt,6)),numel(CCq_g)));
                    fitgmdist_Start(histdata>LeftLim & histdata<RightLim) = pkCt;
                end

                clearvars pks locs w p wx pkParam_Temp pkCt LeftLim RightLim FlagID MinPeakProminenceLim;

                % = Test range of # Gaussian components
                NumGMMComp_Test = 1:NumGMMComp_Est;

                GMM_Cell = cell(1,numel(NumGMMComp_Test));
                GMM_BIC = zeros(1,numel(NumGMMComp_Test)); 
                for gcCt = 1:numel(NumGMMComp_Test) 
                    if isequal(NumGMMComp_Test(gcCt),NumGMMComp_Est)
                        GMM_Cell{1,gcCt} = fitgmdist(histdata,NumGMMComp_Test(gcCt),'CovarianceType','diagonal','RegularizationValue',1E-6,'Start',fitgmdist_Start','Options',GMM_options);
                    else 
                        GMM_Cell{1,gcCt} = fitgmdist(histdata,NumGMMComp_Test(gcCt),'CovarianceType','diagonal','RegularizationValue',1E-6,'Options',GMM_options);
                    end
                    GMM_BIC(1,gcCt)= GMM_Cell{1,gcCt}.BIC;
                end
                
                [minBIC,minBIC_Idx] = min(GMM_BIC); % Smallest BIC value = best
                NumGMMComp = NumGMMComp_Test(1,minBIC_Idx);

                % = Load Best GMM variables
                GMM_mu = GMM_Cell{1,minBIC_Idx}.mu; 
                GMM_sigma = sqrt(permute(GMM_Cell{1,minBIC_Idx}.Sigma,[3 2 1])); 
                GMM_weight = GMM_Cell{1,minBIC_Idx}.ComponentProportion';
                GMM_InGate = (GMM_mu>(tq_g(1)-tq_step)) & (GMM_mu<(tq_g(end)+tq_step)); 

                GMM_mu(~GMM_InGate,:) = []; 
                GMM_sigma(~GMM_InGate,:) = [];
                GMM_weight(~GMM_InGate,:) = [];
                GMM_InGate(~GMM_InGate,:) = [];

                clearvars fitgmdist_Start GMM_options NumGMMComp_Est NumGMMComp_Test GMM_Cell GMM_AIC GMM_BIC minAIC minAIC_Idx minBIC minBIC_Idx NumGMMComp GMM_InGate;

                % === Determine Gaussian component's peak max. height and area
                GMM_pdfMax = zeros(size(GMM_mu),'single'); 
                GMM_pdfArea = zeros(size(GMM_mu),'single'); 
                GMM_pdf = zeros(height(GMM_mu),width(tq_g),'single'); 
                GMM_AllpdfSum = zeros(size(tq_g),'single'); 

                for gcCt = 1:height(GMM_mu) 

                    GMM_pdf(gcCt,:) = normpdf(tq_g,GMM_mu(gcCt),GMM_sigma(gcCt))*GMM_weight(gcCt); 
                    GMM_pdfMax(gcCt,1) = max(GMM_pdf(gcCt,:),[],'all');
                    GMM_AllpdfSum = GMM_AllpdfSum+GMM_pdf(gcCt,:);

                    GMM_cdf = normcdf(tq_g,GMM_mu(gcCt), GMM_sigma(gcCt))*GMM_weight(gcCt); 
                    GMM_pdfArea(gcCt,1) = GMM_cdf(end)-GMM_cdf(1);

                end

                Fig_GMMfit_gCell{gCt,1} = [];        
 
                close all;

                clearvars GMM_pdf GMM_cdf GMM_AllpdfSum cmap gcCt gcCt FigHandle_fr;
                clearvars histdata;
                clearvars IP7_CellCt_GMM; 

                % ===  Load Gaussian component info into GMMParam
                GMMParam_g = ...
                    horzcat(repelem(gCt,numel(GMM_mu))',GMM_mu,GMM_sigma,GMM_weight,GMM_pdfMax,GMM_pdfArea,repelem(Area_g,numel(GMM_mu))'); 
                GMMParam_g = horzcat(GMMParam_g,...
                    GMMParam_g(:,4).*GMMParam_g(:,7)./(1.4826*MDq_full(round(GMMParam_g(:,2))))'); 
                GMMParam_IP7((CellCt_GMM+1):(CellCt_GMM+height(GMMParam_g)),:) = GMMParam_g;

                CellCt_GMM = CellCt_GMM+height(GMMParam_g);

                clearvars tq_g CCq_g MDq_g GMMParam_g GMM_mu GMM_sigma GMM_weight GMM_pdfMax GMM_pdfArea Area_g;

            end
            
            if ~isempty(Fig_GMMfit_gCell)
                Fig_GMMfit_Cell_Out{1,BlkCt} = imtile(Fig_GMMfit_gCell,'GridSize',[NaN 10]); 
            else
                Fig_GMMfit_Cell_Out{1,BlkCt} = []; 
            end

            clearvars Fig_GMMfit_gCell;

            % === Clean up GMMParam
            GMMParam_IP7(GMMParam_IP7(:,1)<0,:) = []; 

            clearvars CellCt_GMM gCt GMM_pdf gcCt CropGate_CCq_full_Plot;


            %% GMM Gaussian components clean up 1: by pdf_thresholding

            clearvars FlagID*;
            
            if ~isempty(GMMParam_IP7)

                % = Tally # Gaussian component found in each CropGate
                Freq_g = tabulate(GMMParam_IP7(:,1));
                Freq_g = Freq_g(:,2); 

                % = Calculate max pdfMax and max Gaussian component ZScore by individual crop-gate
                GMM_GateMaxpdfMax= splitapply(@(x) max(x(:,5),[],'all'),GMMParam_IP7,GMMParam_IP7(:,1));  
                GMM_GateMaxpdfZScore= splitapply(@(x) max(x(:,8),[],'all'),GMMParam_IP7,GMMParam_IP7(:,1));  

                % = Clean up
                % By Each Crop-Gate:
                FlagID_GatepdfMax = zeros(height(GMMParam_IP7),1,'logical'); 
                FlagID_GatepdfZScore = GMMParam_IP7(:,8)<(7.5/18)^2*repelem(GMM_GateMaxpdfZScore,Freq_g); 

                FlagID = FlagID_GatepdfMax | FlagID_GatepdfZScore;

            else 

                FlagID = [];

            end

            GMMParam2_IP7 = GMMParam_IP7;
            GMMParam2_IP7(FlagID,:) = [];

            clearvars Freq_g GMM_GateMaxpdfMax GMM_GateMaxpdfArea GMM_GateMaxpdfZScore GMM_MaxpdfMax GMM_MaxpdfArea;
            clearvars FlagID* gcCt cmap;


            %%  Rebuild CellCt Intensity Profile using pdf_thresholding_cleaned GMM Gaussian Components

            IP8_CellCt = zeros(size(tq_full),'single');
            GMM_pdf = zeros(height(GMMParam2_IP7),width(tq_full),'single');

            for gcCt = 1:height(GMMParam2_IP7) 
                GMM_pdf(gcCt,:) = normpdf(tq_full,GMMParam2_IP7(gcCt,2),GMMParam2_IP7(gcCt,3)).*GMMParam2_IP7(gcCt,4).*GMMParam2_IP7(gcCt,7);
                IP8_CellCt = IP8_CellCt+GMM_pdf(gcCt,:);
            end
            IP8_CellCt = IP8_CellCt'; 
            IP8_CellCt(IP8_CellCt<1E-3*max(IP8_CellCt,[],'all')) = 0; 

            clearvars gcCt GMM_pdf;


            %% Prepare for GMM Gaussian component clean up 2 (w/ findpeaks()): Prepare Crop-Gate

            CropGate_CellCt_IP8 = logical(IP8_CellCt);
            rp_CropGate_CellCt_IP8 = regionprops(logical(CropGate_CellCt_IP8),IP8_CellCt,'Area','PixelIdxList','PixelValues');


            %% GMM Clean Up 2 by findpeaks(): Remove Gaussian Components singly contributing to low peak prominence

            GMMParam_IP9 = GMMParam2_IP7;

            IP9_CellCt = IP8_CellCt;
            CropGate_CellCt_IP9 = CropGate_CellCt_IP8;
            rp_CropGate_CellCt_IP9 = rp_CropGate_CellCt_IP8;


            %% CellCt for IP9: Same as final, but different time step (small tq_step)

            % FPParam_IP9 will hold findpeaks() information
            % Col 1: peak height (pks)
            % Col 2: peak location (locs) 
            % Col 3: peak width (w)  
            % Col 4: peak prominence (p)
            % Col 5: peak width extent, Left Limit  
            % Col 6: peak width extent, Right Limit  
            % Col 7: peak Modified Z-score
            % Col 8: CropGate Index
            % Col 9: peak location (locs)'s tMidPt 

            FPParam_IP9 = -9999*ones(height(GMMParam_IP9),8,'single'); 
            CellCt_IP9 = 0; 
   
            for gCt = 1:height(rp_CropGate_CellCt_IP9) 

                CropGate_Temp = zeros(size(IP9_CellCt),'logical');
                CropGate_Temp(rp_CropGate_CellCt_IP9(gCt).PixelIdxList) = 1;

                % = findpeaks() 
                [IP9_pks,IP9_locs,IP9_p,IP9_w,IP9_wx] = findpeaks_MOD(IP9_CellCt.*CropGate_Temp,'WidthReference','halfheight');

                MinPeakProminenceLim = 0.05*max(IP9_pks+multScMAD_CC*1.4826*IP7_MAD_CC(round(tq_full(IP9_locs))),[],'all'); 
                [IP9_pks,IP9_locs,IP9_p,IP9_w,IP9_wx] = findpeaks_MOD(IP9_CellCt.*CropGate_Temp,'WidthReference','halfheight','MinPeakProminence',MinPeakProminenceLim);
         
                % = Calculate modified Z-score, include in CellCtPkParam (Col 7)
                FPParam_IP9_g = horzcat(IP9_pks,IP9_locs,IP9_w,IP9_p,IP9_wx,... 
                    multScMAD_CC+(IP9_pks./(1.4826*MDq_full(IP9_locs))'),... 
                    repelem(gCt,height(IP9_locs),1)); 

                FPParam_IP9(CellCt_IP9+1:CellCt_IP9+height(IP9_locs),1:8) = FPParam_IP9_g;

                % = Update cell counter
                CellCt_IP9 = CellCt_IP9+height(IP9_locs);

            end

            FPParam_IP9(FPParam_IP9(:,1)<0,:) = []; 

            FlagPairID = find(diff(FPParam_IP9(:,2))<(NumFr_HR/tq_step)); 
            FlagID = zeros(height(FPParam_IP9),1,'logical');
            for pCt = 1:numel(FlagPairID) % pair count
                pks_ratio_Temp = FPParam_IP9(FlagPairID(pCt),1)/FPParam_IP9(FlagPairID(pCt)+1,1);
                if pks_ratio_Temp<(7.5/18)^2
                    FlagID(FlagPairID(pCt),1) = 1;
                elseif pks_ratio_Temp>1/(7.5/18)^2
                    FlagID(FlagPairID(pCt)+1,1) = 1;
                end
            end
            FPParam_IP9(FlagID,:) = [];

            clearvars MinPeakProminenceLim IP9_pks IP9_locs IP9_p IP9_w IP9_wx CellCtPkParam_IP9_g CellCt_IP9 CropGate_Temp gCt FlagPairID FlagID pks_ratio_Temp pCt;


            %% FINAL CellCt

            % CellCtParam_Final will hold CellCt information
            % Col 1: peak height (pks)
            % Col 2: peak location (locs) (@ tq_step = 0.2)
            % Col 3: peak width (w)  (@ tq_step = 0.2)
            % Col 4: peak prominence (p)
            % Col 5: peak width extent, Left Limit  (@ tq_step = 0.2)
            % Col 6: peak width extent, Right Limit  (@ tq_step = 0.2)
            % Col 7: peak Modified Z-score
            % Col 8: CropGate Index
            % Col 9: peak location (locs)'s tMidPt 

            CellCtParam_Final = FPParam_IP9;

            CellCtParam_Final(:,2) = round(tq_full(CellCtParam_Final(:,2))); 
            CellCtParam_Final(:,3) = tq_step.*CellCtParam_Final(:,3); 
            CellCtParam_Final(:,5) = tq_full(round(CellCtParam_Final(:,5))); 
            CellCtParam_Final(:,6) = tq_full(round(CellCtParam_Final(:,6))); 

            % = Include peak location (locs)'s tMidPt
            tMidPt_HistEdge = horzcat(TimeParam_In.tMidPt_FrStart,TimeParam_In.tMidPt_FrEnd(end));
            [tMidPt_HistCt,~,tMidPt_HistBin] = histcounts(CellCtParam_Final(:,2),tMidPt_HistEdge);
            CellCtParam_Final(:,9) = tMidPt_HistBin;
            
            % = Cell Count in each tMidPt
            tMidPt_Final_CellCt_Cell = mat2cell(tMidPt_HistCt,1,repelem(1,numel(tMidPt_HistCt),1));

            % = Mean Z Score in each tMidPt
            if ~isempty(CellCtParam_Final)
                [grp_tMidPt_HistBin,grpID_tMidPt_HistBin] = findgroups(tMidPt_HistBin);
                grp_ZScoreMean = splitapply(@mean, CellCtParam_Final(:,7),grp_tMidPt_HistBin); % # elements = # Groups
                tMidPt_Final_MeanZScore_Cell = zeros(1,TimeParam_In.Num_tMidPt,1);
                tMidPt_Final_MeanZScore_Cell(grpID_tMidPt_HistBin) = grp_ZScoreMean;
                tMidPt_Final_MeanZScore_Cell = mat2cell(tMidPt_Final_MeanZScore_Cell,1,repelem(1,numel(tMidPt_HistCt),1));
            else
                 tMidPt_Final_MeanZScore_Cell = zeros(1,TimeParam_In.Num_tMidPt,1);
                 tMidPt_Final_MeanZScore_Cell = mat2cell(tMidPt_Final_MeanZScore_Cell,1,repelem(1,numel(tMidPt_HistCt),1));
            end
     
            % = Screen report cell count
            fprintf(strcat('\n===========\n'));
            fprintf(strcat('AUTO CellCt:',32,num2str(height(CellCtParam_Final)),'\n'));
            fprintf(strcat('MANUAL CellCt:',32,num2str(numel(Fr_ManualCellCt_In)),' (No Repeat:',32,num2str(numel(unique(Fr_ManualCellCt_In))),')\n'));
            fprintf(strcat('===========\n'));

            clearvars tMidPt_Hist* grp*;


            %% Plot Check: Final Cell Count

            FigHandle = figure('Position',[100 600 1600 2*200+25],'visible','off');
            tlo = tiledlayout(2,1,'TileSpacing','Compact','Padding','tight'); 

            % = Subplot 1
            ax = nexttile(tlo);
            hold(ax, 'on');
            plot((1:1.0:size(LinePlot_In,1))',LinePlot_In,'Color',[0.0 0 0.50],'LineWidth',0.25);
            plot((1:1.0:size(LinePlot_In,1))',IP7_Bkgd,'Color',[0.50 0 0.50],'LineWidth',1.50);

            plot((1:1.0:size(LinePlot_In,1))',IP7_BkgdRmv_pkAug,'Color',[0.5 0 0.50],'LineWidth',1.00); 
            plot((1:1.0:size(LinePlot_In,1))',multScMAD_CC*1.4826*IP7_MAD_CC,'Color',[0.50 0 0.50],'LineWidth',0.25);

            if ~isempty(Fr_ManualCellCt_In)
                xline(Fr_ManualCellCt_In,":",'Color',[0.25 0.75 0],'LineWidth',1.0);
            end
            text(CellCtParam_Final(:,2),LinePlot_In(CellCtParam_Final(:,2))+0.005,'*','Color',[0.25 0 0.25],'FontSize',20,'HorizontalAlignment','Center'); 
            text(CellCtParam_Final(:,2),IP7_BkgdRmv_pkAug(CellCtParam_Final(:,2))+0.005,'*','Color',[0.25 0 0.25],'FontSize',20,'HorizontalAlignment','Center'); 
            xlim([0 height(LinePlot_In)]);
            xticks(0:100:100*floor(0.01*height(LinePlot_In)));
            set(gca,'TickDir','out');
            if ~isempty(Fr_ManualCellCt_In)
                title({strcat('Intensity Profile (IP) Input, IP_BkgdRmv_pkAug'); ...
                    strcat('GMMCellCt:',32,num2str(height(CellCtParam_Final)),'; ',32,...
                    'ManualCellCt:',32,num2str(numel(Fr_ManualCellCt_In)),'; ManualCellCt (no repeat):',32,num2str(numel(unique(Fr_ManualCellCt_In))))},...
                    'FontSize',8,'FontWeight','normal','Interpreter','none');
            else
                title({strcat('Intensity Profile Input, IP_BkgdRmv_pkAug'); ...
                    strcat('GMMCellCt:',32,num2str(height(CellCtParam_Final)))},...
                    'FontSize',8,'FontWeight','normal','Interpreter','none');
            end

            % = Subplot 2
            ax2 = nexttile(tlo);
            hold(ax2, 'on');
            plot(tq_full,IP9_CellCt,'Color',[0.5 0 0.5],'LineWidth',1.00); 
            if ~isempty(Fr_ManualCellCt_In)
                xline(Fr_ManualCellCt_In,":",'Color',[0.25 0.75 0],'LineWidth',1.0);
            end
            text(CellCtParam_Final(:,2),CellCtParam_Final(:,1)+0.0025,'*','Color',[0.25 0 0.25],'FontSize',20,'HorizontalAlignment','Center');
            if ~isempty(Fr_ManualCellCt_In)
                title({strcat('Intensity Profile Input, IP9_CellCt (tq_step=',num2str(tq_step),'fr)'); ...
                    strcat('Peak Z-Score 1 / 25 / 50 / 75 / 99 prctile:',32,strjoin(arrayfun(@num2str,prctile(CellCtParam_Final(:,7),[1 25 50 75 99],'all'),'UniformOutput',false)))},...
                    'FontSize',8,'FontWeight','normal','Interpreter','none');
            else
                title({strcat('Intensity Profile Input, IP9_BkgdRmv_pkAug'); ...
                    strcat('Peak Z-Score 1 / 25 / 50 / 75 / 99 prctile:',32,strjoin(arrayfun(@num2str,prctile(CellCtParam_Final(:,7),[1 25 50 75 99],'all'),'UniformOutput',false)))},...
                    'FontSize',8,'FontWeight','normal','Interpreter','none');
            end
            xlim([0 height(LinePlot_In)]);
            xticks(0:100:100*floor(0.01*height(LinePlot_In)));
            set(gca,'TickDir','out');
           
            Fig_CellCtGMM_Cell_Out{1,BlkCt} = export_fig('-a1','-r150');

            clearvars ax* tlo FigHandle FigHandle_fr;


            %% Create GateXY for each WBC peak in Final CellCt
 
            % i. Locate peaks on IP7_BkgdRmv_pkAug
            [pks,locs,w,p,wx] = findpeaks_MOD(IP7_BkgdRmv_pkAug,'WidthReference','halfprom');
            pkParam_Temp = horzcat(pks,locs,w,p,wx,pks./(1.4826*IP7_MAD_CC(locs))); % Last column (Col7) = Z-Score
            FlagID = pkParam_Temp(:,7)<multScMAD_CC;
            pkParam_Temp(FlagID,:) = [];

            % ii. For each peak in CellCtParam_Final, locate closest peak in
            % pkParam_Temp
            tMidPt_HistEdge = horzcat(TimeParam_In.tMidPt_FrStart,TimeParam_In.tMidPt_FrEnd(end));

            for pkCt = 1:height(CellCtParam_Final)

                [~,~,tMidPt_HistBin] = histcounts(CellCtParam_Final(pkCt,2),tMidPt_HistEdge);
                WinPkWidthEst = tMidPtSkelAv_WinPkWidthEst_Cell{1,tMidPt_HistBin};
               
                [pD,pD_Idx] = pdist2(pkParam_Temp(:,2),CellCtParam_Final(pkCt,2),'euclidean','Smallest',5);
                pD_Idx(pD>WinPkWidthEst) = [];
                pD(pD>WinPkWidthEst) = [];

                GateXY_Temp = zeros(size(LinePlot_In),'logical');
                for iCt = 1:height(pD_Idx)

                    LeftLim = max(1,round(pkParam_Temp(pD_Idx(iCt),2)-2.0/1.1774*(pkParam_Temp(pD_Idx(iCt),2)-pkParam_Temp(pD_Idx(iCt),5))));
                    RightLim = min(round(pkParam_Temp(pD_Idx(iCt),2)+2.0/1.1774*(pkParam_Temp(pD_Idx(iCt),6)-pkParam_Temp(pD_Idx(iCt),2))),height(LinePlot_In));

                    GateXY_Val_Temp = zeros(size(IP7_BkgdRmv_pkAug),'single');
                    GateXY_Val_Temp(LeftLim:RightLim) = IP7_BkgdRmv_pkAug(LeftLim:RightLim);
                    GateXY_Val_Temp(GateXY_Val_Temp<2.0*1.4826*IP7_MAD_CC) = 0;
                    LeftLim = find(GateXY_Val_Temp,1,'first');
                    RightLim = find(GateXY_Val_Temp,1,'last');

                    GateXY_Temp(LeftLim:RightLim) = 1;

                end

                rp = regionprops(GateXY_Temp,IP7_BkgdRmv_pkAug,'Centroid','PixelIdxList','PixelValues','MaxIntensity');

                if pkCt<2
                    rpGateXY = rp;
                else
                    rpGateXY = [rpGateXY; rp];
                end

            end

            clearvars pks locs w p wx pkParam_Temp tMidPt_Hist* FlagID GateXY_Val_Temp GateXY_Temp rp LeftLim RightLim tMidPt_Idx WinPkWidthEst pD pD_Idx pkCt iCt;


            %% Prepare Output

            CellCt_AllBlkorWin(1,BlkCt) = height(CellCtParam_Final);
            CellCt_ExpEr(1,BlkCt) = (TimeParam_In.tMidPt_FrEnd(end)*(1-normcdf(multScMAD_CellCt_In,0,1)))/CellCt_AllBlkorWin(1,BlkCt); 

            CC_pkLoc_AllBlkorWin_Cell{1,BlkCt} = CellCtParam_Final(:,2);
            CC_pkZScore_AllBlkorWin_Cell{1,BlkCt} = CellCtParam_Final(:,7);
            CC_pktMidPt_AllBlkorWin_Cell{1,BlkCt} = CellCtParam_Final(:,9); 
            CC_pkParam_AllBlkorWin_Cell{1,BlkCt} = CellCtParam_Final;

            CC_tMidPt_Ct_AllBlkorWin_Cell{1,BlkCt} = tMidPt_Final_CellCt_Cell; 
            CC_tMidPt_MeanZScore_AllBlkorWin_Cell{1,BlkCt} = tMidPt_Final_MeanZScore_Cell; 
            
            if exist('rpGateXY')
                CC_rpGateXY_AllBlkorWin_Cell{1,BlkCt} = rpGateXY;
            else
                CC_rpGateXY_AllBlkorWin_Cell{1,BlkCt} = single([]);
            end

            CC_IPBkgdRmvPkAg_AllBlkorWin(:,BlkCt) = IP7_BkgdRmv_pkAug;
            CC_IPCC_AllBlkorWin(:,BlkCt) = IP7_CellCt;
            CC_IPMADCC_AllBlkorWin(:,BlkCt) =  IP7_MAD_CC;

        elseif SkipGMM(rCt,BlkCt)

            %% Prepare Output

            Fig_CellCtGMM_Cell_Out{1,BlkCt} = [];
            
            CellCt_AllBlkorWin(1,BlkCt) = height(FPParam_IP7);
            CellCt_ExpEr(1,BlkCt) = (TimeParam_In.tMidPt_FrEnd(end)*(1-normcdf(multScMAD_CellCt_In,0,1)))/CellCt_AllBlkorWin(1,BlkCt);
    
            CC_pkLoc_AllBlkorWin_Cell{1,BlkCt} = FPParam_IP7(:,2);
            CC_pkZScore_AllBlkorWin_Cell{1,BlkCt} = FPParam_IP7(:,7);
            CC_pktMidPt_AllBlkorWin_Cell{1,BlkCt} = FPParam_IP7(:,9); 
            CC_pkParam_AllBlkorWin_Cell{1,BlkCt} = FPParam_IP7;

            CC_tMidPt_Ct_AllBlkorWin_Cell{1,BlkCt} = tMidPt_IP7_CellCt_Cell;
            CC_tMidPt_MeanZScore_AllBlkorWin_Cell{1,BlkCt} = tMidPt_IP7_MeanZScore_Cell;

            CC_rpGateXY_AllBlkorWin_Cell{1,BlkCt} = [];

            CC_IPBkgdRmvPkAg_AllBlkorWin(:,BlkCt) = IP7_BkgdRmv_pkAug;
            CC_IPCC_AllBlkorWin(:,BlkCt) = IP7_CellCt;
            CC_IPMADCC_AllBlkorWin(:,BlkCt) =  IP7_MAD_CC;

        end


    end


    %% Load Output Structure

    CellCtGMM_Struc_Out.CellCt = CellCt_AllBlkorWin;
    CellCtGMM_Struc_Out.Method = repmat("GMM",1,width(CellCt_AllBlkorWin));
    CellCtGMM_Struc_Out.Optimized = ~SkipGMM(rCt,:); 

    CellCtGMM_Struc_Out.CellCtExpEr = CellCt_ExpEr;

    CellCtGMM_Struc_Out.pkLoc = CC_pkLoc_AllBlkorWin_Cell;
    CellCtGMM_Struc_Out.pkZScore = CC_pkZScore_AllBlkorWin_Cell;
    CellCtGMM_Struc_Out.pktMidPt = CC_pktMidPt_AllBlkorWin_Cell;
    CellCtGMM_Struc_Out.pkParam = CC_pkParam_AllBlkorWin_Cell; 
    
    CellCtGMM_Struc_Out.tMidPt_Ct_Cell = CC_tMidPt_Ct_AllBlkorWin_Cell; 
    CellCtGMM_Struc_Out.tMidPt_MeanZScore_Cell = CC_tMidPt_MeanZScore_AllBlkorWin_Cell;

    CellCtGMM_Struc_Out.rpGateXY = CC_rpGateXY_AllBlkorWin_Cell;

    CellCtGMM_Struc_Out.IPCC = CC_IPCC_AllBlkorWin; 
    CellCtGMM_Struc_Out.IPMADCC = CC_IPMADCC_AllBlkorWin; 
    CellCtGMM_Struc_Out.IPBkgdRmvPkAg = CC_IPBkgdRmvPkAg_AllBlkorWin; 

    clearvars CellCt_AllBlkorWin CC*_AllBlkorWin*;


    %% Define Brightness ranking used for deciding which BlkorWin has optimal WBC peak detection

    % = Determine Z-Score cutoff 
    cdf_zs = [-5:0.01:5]; % range of Z-Score to draw cumulative probability function
    cdf_cutoff = (height(LinePlot_BkgdRmv_AllBlkorWin_In)-0.5)/height(LinePlot_BkgdRmv_AllBlkorWin_In);
    cdf_cutoff_zs = cdf_zs(find(cdf('Normal',cdf_zs,0,1)>cdf_cutoff,1,"first"));

    % = Calculate mean Z-score of the top 50% of all Z-Scores > cdf_cutoff_zs
    CellCtGMM_Struc_Out.BrightRankVal = ...
        cellfun(@(x) mean(maxk(x(x>cdf_cutoff_zs),ceil(0.5*numel(x(x>cdf_cutoff_zs))))),CellCtGMM_Struc_Out.pkZScore);

    % = Rank ZScoreHiMean: 1 is the best Blk to use based on WBC peak detection strength
    [~,CellCtGMM_Struc_Out.BrightRank] = sort(CellCtGMM_Struc_Out.BrightRankVal,'descend');

    clearvars cdf_zs cdf_cutoff cdf_cutoff_zs

    if (CtPkGMM_RunOption_In==3) & (rCt<NumRun)
        SkipGMM(NumRun,find(CellCtGMM_Struc_Out.BrightRank==1)) = 0;
    end

end

end
