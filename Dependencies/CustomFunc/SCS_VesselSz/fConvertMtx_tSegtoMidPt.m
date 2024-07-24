function [tMidPt_Val_1i_Mtx_Out,tMidPt_Val_Alli_Mtx_Out] = fConvertMtx_tSegtoMidPt(tSeg_Val_Mtx_In,TimeParam_In)

Numi_In = size(tSeg_Val_Mtx_In,1); % # SkelIdx, SkelPxCt, BlkCt etc
NumSeg_In = TimeParam_In.Num_tSeg; 
NumMidPt_In = TimeParam_In.Num_tMidPt;

% = Sort tFr into tMidPt bins
tMidPt_HistEdge = horzcat(TimeParam_In.tMidPt_FrStart,TimeParam_In.tMidPt_FrEnd(end));
[~,~,tFr_tMidPtBin] = histcounts(1:1:TimeParam_In.tSeg_FrEnd(end),tMidPt_HistEdge);


%% Check if tSeg_Val_Mtx_In has same # Col as NumSeg_In; Match if not

if (width(tSeg_Val_Mtx_In)==1) & (width(tSeg_Val_Mtx_In)<NumSeg_In)
    tSeg_Val_Mtx_ColMatch = repmat(tSeg_Val_Mtx_In,1,NumSeg_In);
else
    tSeg_Val_Mtx_ColMatch = tSeg_Val_Mtx_In;
end


%% Option 1: Treat values from each Dim1 (x) separately
% Output cell size = Numi * NumMidPt

tMidPt_Val_1i_Mtx_Out = NaN(Numi_In,NumMidPt_In,'single'); 

for iCt = 1:Numi_In 

        % = Calculate Segment-averaged Val for One Skeleton
        Val_1i = NaN(TimeParam_In.tSeg_FrEnd(end),NumSeg_In); 
        for SegCt = 1:NumSeg_In 
            Val_1i(TimeParam_In.tSeg_FrStart(SegCt):TimeParam_In.tSeg_FrEnd(SegCt),SegCt) = tSeg_Val_Mtx_ColMatch(iCt,SegCt); % size (tFrCt x NumSeg)
        end
        Val_1i = mean(Val_1i,2,'omitnan'); 

        % = Take mean of values in each tMidPt bin
        tMidPt_Val_1i_Mtx_Out(iCt,:) = (splitapply(@mean,Val_1i,tFr_tMidPtBin'))';

end

clearvars Val_1i iCt SegCt;


%% Option 2: Treat values from all Dim1 elements together
% Output cell size = 1 * NumMidPt
 
tMidPt_Val_Alli_Mtx_Out = NaN(1,NumMidPt_In,'single');

% = Calculate Segment-averaged Val for All Skeletons
Val_Alli = NaN(TimeParam_In.tSeg_FrEnd(end),NumSeg_In,Numi_In); 
for iCt = 1:Numi_In
    for SegCt = 1:NumSeg_In 
        Val_Alli(TimeParam_In.tSeg_FrStart(SegCt):TimeParam_In.tSeg_FrEnd(SegCt),SegCt,iCt) = tSeg_Val_Mtx_ColMatch(iCt,SegCt); 
    end
end

Val_Alli = cell2mat(permute(mat2cell(Val_Alli,height(Val_Alli),width(Val_Alli),repelem(1,Numi_In)),[1 3 2])); 
Val_Alli = mean(Val_Alli,2,'omitnan'); 

% = Take mean of values in each tMidPt bin
tMidPt_Val_Alli_Mtx_Out(1,:) = (splitapply(@mean,Val_Alli,tFr_tMidPtBin'))';

clearvars Val_Alli iCt SegCt;

