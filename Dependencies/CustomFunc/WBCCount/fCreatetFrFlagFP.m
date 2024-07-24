function [tFrFlagFP_Out] = fCreatetFrFlagFP(LinePlot_In,TimeParam_In,ThreshFactor_In,multSigma_In)

tFrFlagFP_Out = zeros(size(LinePlot_In),'logical');

for SegCt=1:TimeParam_In.Num_tSeg

    LeftLim = TimeParam_In.tSeg_FrStart(SegCt);
    RightLim = TimeParam_In.tSeg_FrEnd(SegCt);

    [pks,locs,w,p,wx] = findpeaks_MOD(LinePlot_In,'WidthReference','halfheight'); 
    [npks,nlocs,nw,np,nwx] = findpeaks_MOD(-LinePlot_In,'WidthReference','halfheight'); 
    
    % Keep only peak locations within tSeg 
    % i. Positive peaks (WBC peaks)
    pkParam_Temp = horzcat(pks,locs,w,p,wx);
    FlagID = ~ismember(pkParam_Temp(:,2),LeftLim:RightLim);
    pkParam_Temp(FlagID,:) = [];

    % ii. Negative peaks (RBC pileup troughs)
    npkParam_Temp = horzcat(npks,nlocs,nw,np,nwx);
    FlagID = ~ismember(npkParam_Temp(:,2),LeftLim:RightLim);
    npkParam_Temp(FlagID,:) = [];
    
    % Isolate WBC peaks (positive) and RBC pileup troughs (negative)
    FlagID = ~isoutlier(pkParam_Temp(:,1),"median",'ThresholdFactor',ThreshFactor_In(1)); 
    pkParam_Temp(FlagID,:) = [];

    FlagID = ~isoutlier(npkParam_Temp(:,1),"median",'ThresholdFactor',ThreshFactor_In(2)); 
    npkParam_Temp(FlagID,:) = [];
    
    % Create frame flag for IP2_PrelimBkgd calculation
    for pkCt = 1:height(pkParam_Temp)

        LeftLim = max(1,round(pkParam_Temp(pkCt,2)-multSigma_In/1.1774*(pkParam_Temp(pkCt,2)-pkParam_Temp(pkCt,5)))); 
        RightLim = min(round(pkParam_Temp(pkCt,2)+multSigma_In/1.1774*(pkParam_Temp(pkCt,6)-pkParam_Temp(pkCt,2))),height(LinePlot_In));

        tFrFlagFP_Out(LeftLim:RightLim) = 1;

    end

    for pkCt = 1:height(npkParam_Temp)

        LeftLim = max(1,round(npkParam_Temp(pkCt,2)-multSigma_In/1.1774*(npkParam_Temp(pkCt,2)-npkParam_Temp(pkCt,5))));
        RightLim = min(round(npkParam_Temp(pkCt,2)+multSigma_In/1.1774*(npkParam_Temp(pkCt,6)-npkParam_Temp(pkCt,2))),height(LinePlot_In));

        tFrFlagFP_Out(LeftLim:RightLim) = 1;

    end

end