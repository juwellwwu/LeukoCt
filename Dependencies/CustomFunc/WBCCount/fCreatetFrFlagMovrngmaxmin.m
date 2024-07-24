function [tFrFlagMV_Out] = fCreatetFrFlagMovrngmaxmin(LinePlot_In,movk_In,TimeParam_In,ThreshFactor_In)

% = Calculate movmin, movmax, movrng
if iscell(movk_In)
    IP_movmin = zeros(size(LinePlot_In),'single');
    IP_movmax = zeros(size(LinePlot_In),'single');
    for MidPtCt = 1:TimeParam_In.Num_tMidPt
        IP_movmin_Temp = movmin(LinePlot_In,movk_In{1,MidPtCt});
        IP_movmin(TimeParam_In.tMidPt_FrStart(MidPtCt):TimeParam_In.tMidPt_FrEnd(MidPtCt),1) = ...
            IP_movmin_Temp(TimeParam_In.tMidPt_FrStart(MidPtCt):TimeParam_In.tMidPt_FrEnd(MidPtCt),1);
        IP_movmax_Temp = movmax(LinePlot_In,movk_In{1,MidPtCt});
        IP_movmax(TimeParam_In.tMidPt_FrStart(MidPtCt):TimeParam_In.tMidPt_FrEnd(MidPtCt),1) = ...
            IP_movmax_Temp(TimeParam_In.tMidPt_FrStart(MidPtCt):TimeParam_In.tMidPt_FrEnd(MidPtCt),1);
    end
    IP_movrng =  IP_movmax-IP_movmin;
else
    IP_movmin = movmin(LinePlot_In,movk_In);
    IP_movmax = movmax(LinePlot_In,movk_In);
    IP_movrng =  IP_movmax-IP_movmin;
end

% = Determine tFr to flag by isoutlier()
tFrFlagMovrngmaxmin = ...
    isoutlier(IP_movrng,'median','ThresholdFactor',ThreshFactor_In(1)) | isoutlier(IP_movmax,'median','ThresholdFactor',ThreshFactor_In(2)) | isoutlier(IP_movmin,'median','ThresholdFactor',ThreshFactor_In(3));

% = Dilate such that tFr contributing to flagged tFr's are also flagged
tFrFlagMV_Out = zeros(size(tFrFlagMovrngmaxmin),'logical');
if iscell(movk_In)
    for MidPtCt = 1:TimeParam_In.Num_tMidPt
        FlagID_Temp = movmax(tFrFlagMovrngmaxmin,movk_In{1,MidPtCt});
        tFrFlagMV_Out(TimeParam_In.tMidPt_FrStart(MidPtCt):TimeParam_In.tMidPt_FrEnd(MidPtCt),1) = ...
           FlagID_Temp(TimeParam_In.tMidPt_FrStart(MidPtCt):TimeParam_In.tMidPt_FrEnd(MidPtCt),1);
    end
else
    tFrFlagMV_Out = movmax(tFrFlagMovrngmaxmin,movk_In);
end
