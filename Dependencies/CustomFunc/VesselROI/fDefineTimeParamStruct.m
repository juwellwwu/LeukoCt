function [TimeParam_Struct_Output] = fDefineTimeParamStruct(tSeg_FrLength_Input,ZCropParam_Input)

ImgStack_FrCt_Input = ZCropParam_Input.ZEnd;

% = Define start and end frames for time segments
tSeg_FrSpacing = 0.5*tSeg_FrLength_Input;

tSeg_FrStart = (1:tSeg_FrSpacing:ImgStack_FrCt_Input);
tSeg_FrEnd = tSeg_FrStart+tSeg_FrLength_Input-1;

tSeg_FrStart(tSeg_FrEnd>ImgStack_FrCt_Input)=[];
tSeg_FrEnd(tSeg_FrEnd>ImgStack_FrCt_Input)=[];

tSeg_FrMidPt = round(0.5*(tSeg_FrStart+tSeg_FrEnd));

Num_tSeg = numel(tSeg_FrStart);

% = Define start and end frames for midpoint segments
tMidPt_FrStart = (0.5*tSeg_FrSpacing+1):tSeg_FrSpacing:(ImgStack_FrCt_Input-0.5*tSeg_FrSpacing);
tMidPt_FrEnd = tMidPt_FrStart+tSeg_FrSpacing-1;

tMidPt_FrStart(1) = 1; 
tMidPt_FrEnd(end) = ImgStack_FrCt_Input; 

% = Include Length of tMidPt in frames
tFr_HistEdge = horzcat(tMidPt_FrStart,tMidPt_FrEnd(end));
[tMidPt_FrLength,~,~] = histcounts(1:1:tMidPt_FrEnd(end),tFr_HistEdge);

Num_tMidPt = numel(tMidPt_FrStart);

% = Prepare output structure
TimeParam_Struct_Output = ...
    struct('tSeg_FrLength', single(tSeg_FrLength_Input),...
    'tSeg_FrSpacing', single(tSeg_FrSpacing),...
    'tSeg_FrStart', single(tSeg_FrStart),...
    'tSeg_FrEnd', single(tSeg_FrEnd),...
    'tSeg_FrMidPt', single(tSeg_FrMidPt),...
    'Num_tSeg', single(Num_tSeg),... % = NumSeg 
    'tMidPt_FrStart', single(tMidPt_FrStart),...
    'tMidPt_FrEnd', single(tMidPt_FrEnd),...
    'tMidPt_FrLength', single(tMidPt_FrLength),...
    'Num_tMidPt', single(Num_tMidPt)); % = NumMidPt 
