function [oFlag_RBCBkgd_Out,Fig_RBCBkgd_Flag_Cell_Out] = ...
    fCreateRBCBkgdFlagCell(LinePlot_BkgdRmv_AllBlkorWin_In,LinePlot_MAD_AllBlkorWin_In,LinePlot_Bkgd_AllBlkorWin_In,...
    multScMAD_CellCt_In,...
    tFr_WBCIntensityEst_In,...
    GateXYParam_In,TimeParam_In,SCS_Param_In,NotFrFlag_Option_In)

NumBlk_In = width(LinePlot_BkgdRmv_AllBlkorWin_In);

%% Generate DUMMY Flag

tFr_RBCBkgd_Flag = zeros(size(LinePlot_BkgdRmv_AllBlkorWin_In),'logical');
Fig_RBCBkgd_Flag_Cell_Out = cell(1,NumBlk_In);
oFlag_RBCBkgd_Out = cell(1,NumBlk_In);

for BlkCt = 1:NumBlk_In

    tFr_RBCBkgd_Flag(:,BlkCt) = logical(0);
    Fig_RBCBkgd_Flag_Cell_Out{1,BlkCt} = [];
    
    [oFlag_RBCBkgd_Out{1,BlkCt}] = fConvertFlag_tFrtoMidPtAllSk(tFr_RBCBkgd_Flag(:,BlkCt),0.50,TimeParam_In);

end

