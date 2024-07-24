function [oFlag_FlowStb_Out] = ...
    fCreateDUMMYFlowStabilityFlagCell(oOrient_In,TimeParam_In,NumTimeFr_FlowStbFlag_LoLimit,oSCS_In,ScrnReportBlkChoice_In)

Img_DR_FFT_Filtered_Clean_Cell_In = oOrient_In.Img_DR_OBline_Cell;

NumSkel_In =  oSCS_In.NumSkel;
NumSeg_In =  oSCS_In.NumSeg;
NumMidPt_In =  oSCS_In.NumMidPt;
NumBlk_In =  oSCS_In.NumBlk;


%% Create DUMMY output structure

oFlag_FlowStb_Out = struct;

oFlag_FlowStb_Out.tMidPt_1BlkAllSk_Flag_Cell = repmat({zeros(1,1,'logical')},1,NumMidPt_In,NumBlk_In);  
oFlag_FlowStb_Out.tFr_1BlkAllSk_Flag_Cell = repmat({zeros(TimeParam_In.tMidPt_FrEnd(end),1,'logical')},1,1,NumBlk_In); 

oFlag_FlowStb_Out.tMidPt_1BlkAllSk_Flag_Ct = repmat({zeros(1,1,'logical')},1,1,NumBlk_In);  
oFlag_FlowStb_Out.tMidPt_1BlkAllSk_Flag_Pct = repmat({zeros(1,1,'logical')},1,1,NumBlk_In);

oFlag_FlowStb_Out.tFr_1BlkAllSk_Flag_Ct = repmat({zeros(1,1,'logical')},1,1,NumBlk_In); 
oFlag_FlowStb_Out.tFr_1BlkAllSk_Flag_Pct = repmat({zeros(1,1,'logical')},1,1,NumBlk_In);

oFlag_FlowStb_Out.tFrMidPt_1BlkAllSk_Flag_Ct = repmat({zeros(1,1,'logical')},1,1,NumBlk_In); 
oFlag_FlowStb_Out.tFrMidPt_1BlkAllSk_Flag_Pct = repmat({zeros(1,1,'logical')},1,1,NumBlk_In);

oFlag_FlowStb_Out.tMidPt_AllBlkAllSk_Flag_Ct = 0; 
oFlag_FlowStb_Out.tMidPt_AllBlkAllSk_Flag_Pct = 0; 

oFlag_FlowStb_Out.tFr_AllBlkAllSk_Flag_Ct = 0; 
oFlag_FlowStb_Out.tFr_AllBlkAllSk_Flag_Pct = 0;

oFlag_FlowStb_Out.tFrMidPt_AllBlkAllSk_FlagCt = 0; 
oFlag_FlowStb_Out.tFrMidPt_AllBlkAllSk_FlagPct = 0;

oFlag_FlowStb_Out.tFrNonFlag_AllBlkAllSk_MaxContCt = TimeParam_In.tMidPt_FrEnd(end); 
oFlag_FlowStb_Out.tFrNonFlag_AllBlkAllSk_MaxContCt_FrStart = 1;
oFlag_FlowStb_Out.tFrNonFlag_AllBlkAllSk_MaxContCt_FrEnd = TimeParam_In.tMidPt_FrEnd(end);





