function [LinePlot_Bkgd_AllBlkorWin_Out,LinePlot_BkgdRmv_AllBlkorWin_Out,LinePlot_MAD_AllBlkorWin_Out,Fig_Bkgd_BkgdRmv_MAD_Cell_Out] = ...
    fAnalyzeFrIntensityProfile(LinePlot_AllBlkorWin_In, ...
    multScMAD_CellCt_In,... 
    GateXYParam_In,TimeParam_In,XY_PxLength_In,Z_PxLength_In,...
    Velocity_Param_In,...  
    Fr_ManualCellCt_In)  

NumBlk_In = width(LinePlot_AllBlkorWin_In);

multScMAD_FP_pk = 2.65; 
multScMAD_MV_pk = 3.00; 
multScMAD_CC = 3.50; 

tMidPtSkelAv_Velocity_Cell_In = Velocity_Param_In.tMidPt_AllSkOBwMean_Cell; 


%% # Time Frames For Smoothing

% 1) # tFr in 1 sec
% Only used for preliminary smoothing
NumFr_1sec = 1/Z_PxLength_In;

% 2) # tFr based on human heart rate
NumFr_HR = round((60/Z_PxLength_In)/240); % default: round((60/Z_PxLength_In)/240); ** must check all step's accuracy if change **


%%

LinePlot_Bkgd_AllBlkorWin_Out = zeros(size(LinePlot_AllBlkorWin_In),'single');
LinePlot_BkgdRmv_AllBlkorWin_Out = zeros(size(LinePlot_AllBlkorWin_In),'single');
LinePlot_MAD_AllBlkorWin_Out = zeros(size(LinePlot_AllBlkorWin_In),'single');

Fig_Bkgd_BkgdRmv_MAD_Cell_Out = cell(1,NumBlk_In);

for BlkCt = 1:NumBlk_In


    LinePlot_In = LinePlot_AllBlkorWin_In(:,BlkCt);


    %% GateXYParam, from fLocateOptimalGateXYWindow()

    WBCtoGateXYWinAreaRatio_In = GateXYParam_In.WBCtoGateXYWinAreaRatio(:,BlkCt);
    WinLength_px_In = GateXYParam_In.WinLength_px(:,BlkCt);


    %% Expected WBC peak width in tFr

    tMidPtSkelAv_WinPkWidthEst_Cell = cellfun(@(x) WinLength_px_In/x,tMidPtSkelAv_Velocity_Cell_In,'UniformOutput',false);


    %% 1) Preliminary peak loc estimate

    IP1 = LinePlot_In;

    % = Perform very rough background removal
    IP1_RoughBkgd= movmedian(filloutliers(IP1,"center","movmedian",NumFr_1sec),NumFr_HR);
    IP1_RoughBkgdRmv = IP1-IP1_RoughBkgd;

    % = Flag WBC peaks and RBC pileup troughs
    IP1_tFrFlag = ...
        fCreatetFrFlagMovrngmaxmin(IP1_RoughBkgdRmv,cellfun(@(x) round(3*x), tMidPtSkelAv_WinPkWidthEst_Cell,'UniformOutput',false),TimeParam_In,[multScMAD_MV_pk multScMAD_MV_pk multScMAD_MV_pk]);

    clearvars IP1_RoughBkgd IP1_RoughBkgdRmv IP1_tFrFlag_Plot ThreshFactor;


    %% 2) Create preliminary background, subtract from LinePlot_In

    % = Replace rough peak locations from IP1 with NaN, then linear interpolate
    % background in NaN regions
    IP2 = LinePlot_In;
    IP2(IP1_tFrFlag) = NaN; 
    IP2 = fFillMissingLinMethod(IP2);

    % = Smooth to generate prelim background
    sgolay_FrLength = NumFr_HR;
    sgolay_Order = max(3,round(sgolay_FrLength/(max(cell2mat(tMidPtSkelAv_WinPkWidthEst_Cell),[],'all')*5))); 

    IP2_PrelimBkgd = fCreateFrIntensityProfileBkgdSG(IP2,sgolay_Order,sgolay_FrLength,Z_PxLength_In,'n',2);
    IP2_PrelimBkgdRmv = LinePlot_In-IP2_PrelimBkgd;

    clearvars IP_movrng_FrFlagID;
    clearvars sgolay_FrLength sgolay_Order;


    %% 3) Create final background: i. Repeat peak loc estimate, this time using background-removed (flattened) Prelim intensity profile

    IP3 = IP2_PrelimBkgdRmv; 

    % = Create tFr Flag for IP4_MAD Creation
    IP3_FP_tFrFlag = fCreatetFrFlagFP(IP3,TimeParam_In,[multScMAD_FP_pk multScMAD_FP_pk],4.0);
    IP3_MV_tFrFlag = fCreatetFrFlagMovrngmaxmin(IP3,cellfun(@(x) round(3*x), tMidPtSkelAv_WinPkWidthEst_Cell,'UniformOutput',false),TimeParam_In,[multScMAD_MV_pk multScMAD_MV_pk multScMAD_MV_pk]);

    IP3_tFrFlag = IP3_FP_tFrFlag | IP3_MV_tFrFlag;

    clearvars IP3_tFrFlag_Plot ThreshFactor;
    clearvars IP3; 


    %% 4) Create final background: ii. Measure time-dependent "noise" of background-removed (flattened) prelim intensity profile

    IP4 = IP2_PrelimBkgdRmv;
    IP4_MAD = fCalcMovFun(IP4, "movmad", NumFr_HR, 2.00, IP3_tFrFlag);

    sgolay_FrLength = NumFr_HR;
    sgolay_Order = max(3,round(sgolay_FrLength/(max(cell2mat(tMidPtSkelAv_WinPkWidthEst_Cell),[],'all')*5))); 
    IP4_MAD = fCreateFrIntensityProfileBkgdSG(IP4_MAD,sgolay_Order,sgolay_FrLength,Z_PxLength_In,'n',2);

    clearvars IP4 IP3_tFrFlag_Plot; 


    %% Create final background, final background-removed (flattened) Intensity profile, final (time-dependent) MAD
 
    for itCt = 1:3

        if itCt==1
            IP_PrelimBkgd_In = IP2_PrelimBkgd;
            IP_PrelimBkgdRmv_In = IP2_PrelimBkgdRmv;
            IP_PrelimMAD_In = IP4_MAD;
            IP_MV_tFrFlag_In = IP3_MV_tFrFlag;
        elseif itCt>1
            IP_PrelimBkgd_In = IP7_Bkgd;
            IP_PrelimBkgdRmv_In = IP7_BkgdRmv;
            IP_PrelimMAD_In = IP7_MAD;
            IP_MV_tFrFlag_In = IP7_MV_tFrFlag;
        end

        %% 5a) Create final background: iiia. Use findpeaks() to locate WBC peak loc in "flattened" IP2_PrelimBkgdRmv

        % i. WBC peaks bright: peaks in upward direction

        IP5 = IP_PrelimBkgdRmv_In; 

        % =  IP5_pkLoc_Gate specifies where locs will be allowed
        IP5_IP6_pk_multScMAD = multScMAD_FP_pk;
        IP5_pkLoc_Gate = IP5-(IP5_IP6_pk_multScMAD*1.4826*IP_PrelimMAD_In);
        IP5_pkLoc_Gate(IP5_pkLoc_Gate<0) = 0;
        IP5_pkLoc_Gate = logical(IP5_pkLoc_Gate);

        [IP5_pks,IP5_locs,IP5_w,IP5_p,IP5_wx] = findpeaks_MOD(IP5,'WidthReference','halfheight');

        IP5_pkFlagID = (~ismember(IP5_locs,find(IP5_pkLoc_Gate))) & (~ismember(IP5_locs,find(IP_MV_tFrFlag_In)));

        IP5_pkParam = horzcat(IP5_pks,IP5_locs,IP5_p,IP5_w,IP5_wx);
        IP5_pkParam(IP5_pkFlagID,:)  = [];

        clearvars IP5 IP5_pkLoc_Gate IP5_pkFlagID IP5_pks IP5_locs IP5_p IP5_w IP5_wx; 


        %% 5b) Create final background: iiib. Use findpeaks() to locate RBC pileup trough in "flattened" IP2_PrelimBkgdRmv

        % ii. RBC pileups dark: peaks in downward direction

        IP6 = (-IP_PrelimBkgdRmv_In); 

        % =  IP6_pkLoc_Gate specifies where locs will be allowed
        IP6_pkLoc_Gate = IP6-(IP5_IP6_pk_multScMAD*1.4826*IP_PrelimMAD_In);
        IP6_pkLoc_Gate(IP6_pkLoc_Gate<0) = 0;
        IP6_pkLoc_Gate = logical(IP6_pkLoc_Gate);

        [IP6_pks,IP6_locs,IP6_w,IP6_p,IP6_wx] = findpeaks_MOD(IP6,'WidthReference','halfheight');

        IP6_pkFlagID = (~ismember(IP6_locs,find(IP6_pkLoc_Gate))) & (~ismember(IP6_locs,find(IP_MV_tFrFlag_In)));

        IP6_pkParam = horzcat(IP6_pks,IP6_locs,IP6_p,IP6_w,IP6_wx); 
        IP6_pkParam(IP6_pkFlagID,:)  = [];

        clearvars IP6 IP6_pkLoc_Gate IP6_pkFlagID IP6_pks IP6_locs IP6_p IP6_w IP6_wx; 


        %% 6) Create final background: iv. Build WBC Peak & RBC pileup trough location Flag, for final background and MAD calculation

        multSigma = [2.0 3.5 5.0];

        % i) Flag with WBC peaks only
        tFr_pk_Up_Flag = zeros(height(LinePlot_In),numel(multSigma),'logical');
        for mCt = 1:numel(multSigma) 
            for pkCt=1:height(IP5_pkParam)
                LeftLim = max(1,floor(IP5_pkParam(pkCt,2)-multSigma(mCt)/1.1774*(IP5_pkParam(pkCt,2)-IP5_pkParam(pkCt,5))));
                RightLim = min(ceil(IP5_pkParam(pkCt,2)+multSigma(mCt)/1.1774*(IP5_pkParam(pkCt,6)-IP5_pkParam(pkCt,2))),height(LinePlot_In));
                tFr_pk_Up_Flag(LeftLim:RightLim,mCt) = 1;
            end
        end

        % ii) Flag with RBC pileup troughs only
        tFr_pk_Dn_Flag = zeros(height(LinePlot_In),numel(multSigma),'logical');
        for mCt = 1:numel(multSigma) 
            for pkCt=1:height(IP6_pkParam)
                LeftLim = max(1,floor(IP6_pkParam(pkCt,2)-multSigma(mCt)/1.1774*(IP6_pkParam(pkCt,2)-IP6_pkParam(pkCt,5))));
                RightLim = min(ceil(IP6_pkParam(pkCt,2)+multSigma(mCt)/1.1774*(IP6_pkParam(pkCt,6)-IP6_pkParam(pkCt,2))),height(LinePlot_In));
                tFr_pk_Dn_Flag(LeftLim:RightLim,mCt) = 1;
            end
        end

        % iii) Flag with both WBC peaks and RBC pileup troughs (for MAD calculation)
        tFr_pk_Flag = zeros(height(LinePlot_In),numel(multSigma),'logical');
        for mCt = 1:numel(multSigma) 
            tFr_pk_Flag(:,mCt) = max(tFr_pk_Up_Flag(:,mCt),tFr_pk_Dn_Flag(:,mCt));
        end

        clearvars pkCt mCt LeftLim RightLim multFWHM tFr_pk_*_Flag_Plot IP_*_tFrFlag_In_Plot tFr_pk_Flag_Plot;


        %%  7) Create final background: v.) Create final (time-dependent) background of LinePlot_In
 
        IP7 = LinePlot_In;

        IP7_Bkgd = IP7;
        IP7_Bkgd(tFr_pk_Flag(:,multSigma==2)|IP_MV_tFrFlag_In) = NaN; 

        IP7_Bkgd_Flg = IP7;
        IP7_Bkgd_Flg(IP7_Bkgd_Flg>(IP_PrelimBkgd_In+IP5_IP6_pk_multScMAD*1.4826*IP_PrelimMAD_In) | IP7_Bkgd_Flg<(IP_PrelimBkgd_In-IP5_IP6_pk_multScMAD*1.4826*IP_PrelimMAD_In)) = NaN;

        IP7_Bkgd(tFr_pk_Flag(:,multSigma==2)|IP_MV_tFrFlag_In) = IP7_Bkgd_Flg(tFr_pk_Flag(:,multSigma==2)|IP_MV_tFrFlag_In);

        IP7_Bkgd = fFillMissingLinMethod(IP7_Bkgd); 

        sgolay_FrLength = NumFr_HR;
        sgolay_Order = max(3,round(sgolay_FrLength/(max(cell2mat(tMidPtSkelAv_WinPkWidthEst_Cell),[],'all')*5))); 
        IP7_Bkgd = fCreateFrIntensityProfileBkgdSG(IP7_Bkgd,sgolay_Order,sgolay_FrLength,Z_PxLength_In,'n',2);

        clearvars IP7_Bkgd_Flg;
        clearvars pkCt LeftLim RightLim sgolay_FrLength sgolay_Order;


        %%  8) Create final background-removed (flattened) Intensity Profile

        IP7_BkgdRmv = IP7-IP7_Bkgd;

        clearvars IP7; 


        %% 9) Locate large range areas in final background-removed (flattened) Intensity Profile for MAD Calculation

        IP7_MV_tFrFlag = ...
            fCreatetFrFlagMovrngmaxmin(IP7_BkgdRmv,cellfun(@(x) round(3*x), tMidPtSkelAv_WinPkWidthEst_Cell,'UniformOutput',false),TimeParam_In,[multScMAD_MV_pk multScMAD_MV_pk multScMAD_MV_pk]);

        clearvars movk IP7_movmin IP7_movmax IP7_movrng ThreshFactor;


        %%  10) Create final MAD of Intensity Profile for Cell Counting

        tFr_SG_Flag = tFr_pk_Flag(:,multSigma==5)|IP7_MV_tFrFlag;
        IP7_MAD = fCalcMovFun(IP7_BkgdRmv, "movmad",NumFr_HR, 2.00, tFr_SG_Flag);

        clearvars sgolay_FrLength sgolay_Order;


        %% Plot Check

        FigHandle = figure('Position',[100 600 1600 1*200+25],'visible','off'); hold on;

        plot((1:1.0:size(LinePlot_In,1))',LinePlot_In,'Color',[0.0 0 0.50],'LineWidth',0.25);
        plot((1:1.0:size(LinePlot_In,1))',IP7_Bkgd,'Color',[0.0 0 0.50],'LineWidth',1.00);
        plot((1:1.0:size(LinePlot_In,1))',IP7_BkgdRmv,'Color',[0.5 0 0.50],'LineWidth',1.00);
        plot((1:1.0:size(LinePlot_In,1))',multScMAD_CC*1.4826*IP7_MAD,'Color',[0.50 0 0.50],'LineWidth',0.25);
        plot((1:1.0:size(LinePlot_In,1))',-multScMAD_CC*1.4826*IP7_MAD,'Color',[0.50 0 0.50],'LineWidth',0.25);
        if ~isempty(Fr_ManualCellCt_In)
            xline(Fr_ManualCellCt_In,":",'Color',[0.25 0.75 0],'LineWidth',1.0);
        end
        xlim([0 height(LinePlot_In)]);
        xticks(0:100:100*floor(0.01*height(LinePlot_In)));
        set(gca,'TickDir','out');
        if ~isempty(Fr_ManualCellCt_In)
            title({strcat('Blue: LinePlot_In (BlkorWin=',32,num2str(BlkCt),'); Purple: LinePlot_Bkgd+multScMAD_CC*1.4826*LinePlot_MAD'); ...
                strcat('ManualCellCt:',32,num2str(numel(Fr_ManualCellCt_In)),'; ManualCellCt (no repeat):',32,num2str(numel(unique(Fr_ManualCellCt_In))))},...
                'FontSize',8,'FontWeight','normal','Interpreter','none');
        else
            title({strcat('Blue: LinePlot_In (BlkorWin=',32,num2str(BlkCt),'; Purple: LinePlot_Bkgd+multScMAD_CC*1.4826*LinePlot_MAD')},...
                'FontSize',8,'FontWeight','normal','Interpreter','none');
        end

        Fig_Bkgd_BkgdRmv_MAD_Cell_Out{1,BlkCt} = export_fig('-a1','-r150');

        clearvars FigHandle_fr;


        %% Load in output

        LinePlot_Bkgd_AllBlkorWin_Out(:,BlkCt) = IP7_Bkgd;
        LinePlot_BkgdRmv_AllBlkorWin_Out(:,BlkCt) = IP7_BkgdRmv;
        LinePlot_MAD_AllBlkorWin_Out(:,BlkCt) = IP7_MAD;

    end

    clearvars IP_PrelimBkgd_In IP_PrelimBkgdRmv_In IP_PrelimMAD_In IP_*_tFrFlag_In tCt tFr_pk_*Flag tFr_SG_Flag multSigma itCt;


end

close all;


%% Clear obsolete variables

clearvars *multScMAD*;
clearvars tMidPtSkelAv_EstWBCPkWidth_Cell;
clearvars CropGate* rp_CropGate*;
clearvars GMMOverlap*;
clearvars t_full tq_full MDq_full CCq_full rp_CCq_full;
clearvars NumFr_1sec NumFr_HR;
clearvars LinePlot_In WBCtoGateXYWinAreaRatio_In WinLength_px_In tMidPtSkelAv_Velocity_Cell_In XY_PxLength_In Z_PxLength_In TimeParam_In NumBlk_In;

end


%% =========
%% ========= IN SCRIPT FUNCTIONS

%%

function [LinePlot_MovFun_Out] = fCalcMovFun(LinePlot_In, MovFunOption, MovFun_k_In, MovFun_k_Tol_In, Px_Flag_In)

if isrow(LinePlot_In)
    LinePlot = LinePlot_In';
else
    LinePlot = LinePlot_In;
end

% = Determine the indices for each pixel for movmedian
pD_In = (1:1:height(LinePlot))';
pD_In(Px_Flag_In) = height(LinePlot).*10; 

[pD,pD_Idx] = pdist2(pD_In,(1:1:height(LinePlot))','euclidean','Smallest',MovFun_k_In);
pD_Idx(pD>MovFun_k_Tol_In*0.5*MovFun_k_In) = NaN; 
pD(pD>MovFun_k_Tol_In*0.5*MovFun_k_In) = NaN;

tFr_MovCalc_Cell =  mat2cell(pD_Idx,height(pD_Idx),repelem(1,width(pD_Idx)));
tFr_MovCalc_Cell = cellfun(@(x) x(~isnan(x)),tFr_MovCalc_Cell,'UniformOutput',false);

if MovFunOption == 'movmad'
    LinePlot_MovFun_Out = cell2mat(cellfun(@(x) mad(LinePlot(x)),tFr_MovCalc_Cell,'UniformOutput',false)); 
elseif MovFunOption == 'movmedian'
    LinePlot_MovFun_Out = cell2mat(cellfun(@(x) median(LinePlot(x)),tFr_MovCalc_Cell,'UniformOutput',false)); 
end

LinePlot_MovFun_Out(cellfun(@(x) height(x),tFr_MovCalc_Cell,'UniformOutput',true)<max(round(0.50*MovFun_k_In),30)) = NaN;

LinePlot_MovFun_Out  = fFillMissingLinMethod(LinePlot_MovFun_Out);

end


%% 

function [LinePlot_FillMissingLinMethod_Out] = fFillMissingLinMethod(LinePlot_In)

if iscolumn(LinePlot_In)
    LinePlot_FillMissingLinMethod_Out = LinePlot_In';
elseif isrow(LinePlot_In)
    LinePlot_FillMissingLinMethod_Out = LinePlot_In;
end

LinePlot_FillMissingLinMethod_Out = ...  
    horzcat(LinePlot_FillMissingLinMethod_Out(find(~isnan(LinePlot_FillMissingLinMethod_Out),1,'first')),LinePlot_FillMissingLinMethod_Out,LinePlot_FillMissingLinMethod_Out(find(~isnan(LinePlot_FillMissingLinMethod_Out),1,'last')));
LinePlot_FillMissingLinMethod_Out = fillmissing(LinePlot_FillMissingLinMethod_Out,'linear');
LinePlot_FillMissingLinMethod_Out = (LinePlot_FillMissingLinMethod_Out(2:(end-1)))';

end
