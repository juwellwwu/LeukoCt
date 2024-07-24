function [oFlag_VctyEr_Out] = fCreateVelocityStdErFlagCell(oVelocity_In,TimeParam_In,oSCS_In,StdErLimOB_In,StdErLimSK_In,ScrnReportBlkChoice_In)

NumMidPt_In = oSCS_In.NumMidPt;
NumBlk_In = oSCS_In.NumBlk;

% Assign tFr to different tMidPt
tFr_HistEdge = horzcat(TimeParam_In.tMidPt_FrStart,TimeParam_In.tMidPt_FrEnd(end));
[tFr_HistCt,~,tFr_HistBin] = histcounts(1:1:TimeParam_In.tMidPt_FrEnd(end),tFr_HistEdge);

oFlag_VctyEr_Out = struct;


%% Calculate Velocity error: i. Velocity Error from Orient Band Statistics

% Per Vessel Block
tMidPt_Blk_Metric_AllSkOBwStdMean_Cell = cell2mat(oVelocity_In.tMidPt_AllSkOBwStdEr_Cell).*2./cell2mat(oVelocity_In.tMidPt_AllSkOBwMean_Cell); 
tMidPt_Blk_Flag_AllSkOBwStdMean_Cell = tMidPt_Blk_Metric_AllSkOBwStdMean_Cell > StdErLimOB_In;

% If oVelocity.tMidPt_AllSkOBct_Cell = 0 or 1; ie, if tMidPt doesn't
% have orient bands (and orient and velocity values are interpolated), also
% flag the cell. 
tMidPt_Blk_Flag_AllSkOBwStdMean_Cell = tMidPt_Blk_Flag_AllSkOBwStdMean_Cell | (cell2mat(oVelocity_In.tMidPt_AllSkOBct_Cell)<2);

tMidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell = zeros(1,1,NumBlk_In,'single'); 
tMidPt_Blk_Flag_AllSkOBwStdMean_Pct_Cell = zeros(1,1,NumBlk_In,'single'); 

tMidPtwMean_Blk_Metric_AllSkOBwStdMean_Cell = zeros(1,1,NumBlk_In,'single'); 
tMidPtwMedian_Blk_Metric_AllSkOBwStdMean_Cell = zeros(1,1,NumBlk_In,'single'); 

tFr_Blk_Flag_AllSkOBwStdMean_Cell = zeros(TimeParam_In.tMidPt_FrEnd(end),1,NumBlk_In,'logical');
tFr_Blk_Flag_AllSkOBwStdMean_Ct_Cell = zeros(1,1,NumBlk_In,'single'); 
tFr_Blk_Flag_AllSkOBwStdMean_Pct_Cell = zeros(1,1,NumBlk_In,'single');

for BlkCt = 1:NumBlk_In

    tMidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell(1,1,BlkCt) = sum(tMidPt_Blk_Flag_AllSkOBwStdMean_Cell(1,:,BlkCt),2);
    tMidPt_Blk_Flag_AllSkOBwStdMean_Pct_Cell(1,1,BlkCt) = tMidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell(1,1,BlkCt)/NumMidPt_In*100;

    tMidPtwMean_Blk_Metric_AllSkOBwStdMean_Cell(1,1,BlkCt) = ...
        sum(TimeParam_In.tMidPt_FrLength.*tMidPt_Blk_Metric_AllSkOBwStdMean_Cell(1,:,BlkCt),2)/sum(TimeParam_In.tMidPt_FrLength,2);

    % weightedMedian() from FileExchange:
    % https://www.mathworks.com/matlabcentral/fileexchange/23077-weighted-median
    tMidPtwMedian_Blk_Metric_AllSkOBwStdMean_Cell(1,1,BlkCt) = ...
        weightedMedian(tMidPt_Blk_Metric_AllSkOBwStdMean_Cell(1,:,BlkCt),TimeParam_In.tMidPt_FrLength);

    tFr_Flag_Temp = zeros(TimeParam_In.tMidPt_FrEnd(end),1,1,'logical');
    tFr_Flag_Temp(ismember(tFr_HistBin,find(tMidPt_Blk_Flag_AllSkOBwStdMean_Cell(:,:,BlkCt))),1) = 1;
    tFr_Blk_Flag_AllSkOBwStdMean_Cell(:,1,BlkCt) = tFr_Flag_Temp;

    tFr_Blk_Flag_AllSkOBwStdMean_Ct_Cell(1,1,BlkCt) = sum(tFr_HistCt(find(tMidPt_Blk_Flag_AllSkOBwStdMean_Cell(:,:,BlkCt))),'all');
    tFr_Blk_Flag_AllSkOBwStdMean_Pct_Cell(1,1,BlkCt) = tFr_Blk_Flag_AllSkOBwStdMean_Ct_Cell(1,1,BlkCt)/TimeParam_In.tMidPt_FrEnd(end)*100;

end

% All vessel blocks
tMidPt_Flag_AllSkOBwStdMean_Ct = sum(tMidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell,'all'); 
tMidPt_Flag_AllSkOBwStdMean_Pct = tMidPt_Flag_AllSkOBwStdMean_Ct/(NumBlk_In*NumMidPt_In)*100; 

tMidPt_Metric_AllSkOBwStdMean_Mean = mean(tMidPt_Blk_Metric_AllSkOBwStdMean_Cell,'all'); 
tMidPt_Metric_AllSkOBwStdMean_Median = median(tMidPt_Blk_Metric_AllSkOBwStdMean_Cell,'all'); 

tFr_Flag_AllSkOBwStdMean_Ct = sum(tFr_Blk_Flag_AllSkOBwStdMean_Ct_Cell,'all'); 
tFr_Flag_AllSkOBwStdMean_Pct = tFr_Flag_AllSkOBwStdMean_Ct/(NumBlk_In*TimeParam_In.tMidPt_FrEnd(end))*100;

% = Determine longest non-flagged stretch of time (for all vessel blocks)
rp_NonFlag_AllSkOBwStdMean = regionprops(~logical(max(tFr_Blk_Flag_AllSkOBwStdMean_Cell,[],3)),'Area','PixelIdxList');

if ~isempty(rp_NonFlag_AllSkOBwStdMean) 
    [MaxNonFlagtFrCt_AllSkOBwStdMean, MaxNonFlagtFrCt_AllSkOBwStdMean_BlobIdx] = max(cat(1,rp_NonFlag_AllSkOBwStdMean.Area),[],'all');
    MaxNonFlagtFrCt_AllSkOBwStdMean_FrStart = rp_NonFlag_AllSkOBwStdMean(MaxNonFlagtFrCt_AllSkOBwStdMean_BlobIdx).PixelIdxList(1);
    MaxNonFlagtFrCt_AllSkOBwStdMean_FrEnd = rp_NonFlag_AllSkOBwStdMean(MaxNonFlagtFrCt_AllSkOBwStdMean_BlobIdx).PixelIdxList(end);
else 
    MaxNonFlagtFrCt_AllSkOBwStdMean = 0;
    MaxNonFlagtFrCt_AllSkOBwStdMean_FrStart = 0;
    MaxNonFlagtFrCt_AllSkOBwStdMean_FrEnd = 0;
end

% = Convert to cell
tMidPt_Blk_Metric_AllSkOBwStdMean_Cell = mat2cell(tMidPt_Blk_Metric_AllSkOBwStdMean_Cell,1,repelem(1,NumMidPt_In),repelem(1,NumBlk_In));
tMidPt_Blk_Flag_AllSkOBwStdMean_Cell = mat2cell(tMidPt_Blk_Flag_AllSkOBwStdMean_Cell,1,repelem(1,NumMidPt_In),repelem(1,NumBlk_In));
tMidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell = mat2cell(tMidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell,1,1,repelem(1,NumBlk_In));
tMidPt_Blk_Flag_AllSkOBwStdMean_Pct_Cell = mat2cell(tMidPt_Blk_Flag_AllSkOBwStdMean_Pct_Cell,1,1,repelem(1,NumBlk_In));

tMidPtwMean_Blk_Metric_AllSkOBwStdMean_Cell = mat2cell(tMidPtwMean_Blk_Metric_AllSkOBwStdMean_Cell,1,1,repelem(1,NumBlk_In));
tMidPtwMedian_Blk_Metric_AllSkOBwStdMean_Cell = mat2cell(tMidPtwMedian_Blk_Metric_AllSkOBwStdMean_Cell,1,1,repelem(1,NumBlk_In));

tFr_Blk_Flag_AllSkOBwStdMean_Cell = mat2cell(tFr_Blk_Flag_AllSkOBwStdMean_Cell,TimeParam_In.tMidPt_FrEnd(end),1,repelem(1,NumBlk_In));
tFr_Blk_Flag_AllSkOBwStdMean_Ct_Cell = mat2cell(tFr_Blk_Flag_AllSkOBwStdMean_Ct_Cell,1,1,repelem(1,NumBlk_In));
tFr_Blk_Flag_AllSkOBwStdMean_Pct_Cell = mat2cell(tFr_Blk_Flag_AllSkOBwStdMean_Pct_Cell,1,1,repelem(1,NumBlk_In));

tFr_MidPt_Blk_Flag_AllSkOBwStdMean_Cell = tFr_Blk_Flag_AllSkOBwStdMean_Cell;
tFr_MidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell = tFr_Blk_Flag_AllSkOBwStdMean_Ct_Cell;
tFr_MidPt_Blk_Flag_AllSkOBwStdMean_Pct_Cell = tFr_Blk_Flag_AllSkOBwStdMean_Pct_Cell;
tFr_MidPt_Flag_AllSkOBwStdMean_Ct = tFr_Flag_AllSkOBwStdMean_Ct;
tFr_MidPt_Flag_AllSkOBwStdMean_Pct = tFr_Flag_AllSkOBwStdMean_Pct; 

% = Load into struct
oFlag_VctyEr_Out.OB.tMidPt_1BlkAllSk_Flag_Cell = tMidPt_Blk_Flag_AllSkOBwStdMean_Cell;
oFlag_VctyEr_Out.OB.tMidPt_1BlkAllSk_Metric_Cell = tMidPt_Blk_Metric_AllSkOBwStdMean_Cell; 

oFlag_VctyEr_Out.OB.tFr_1BlkAllSk_Flag_Cell = tFr_Blk_Flag_AllSkOBwStdMean_Cell; 

oFlag_VctyEr_Out.OB.tMidPt_1BlkAllSk_Flag_Ct = tMidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell; 
oFlag_VctyEr_Out.OB.tMidPt_1BlkAllSk_Flag_Pct = tMidPt_Blk_Flag_AllSkOBwStdMean_Pct_Cell;

oFlag_VctyEr_Out.OB.tMidPtwMean_1BlkAllSk_Metric_Cell = tMidPtwMean_Blk_Metric_AllSkOBwStdMean_Cell;
oFlag_VctyEr_Out.OB.tMidPtwMedian_1BlkAllSk_Metric_Cell = tMidPtwMedian_Blk_Metric_AllSkOBwStdMean_Cell;

oFlag_VctyEr_Out.OB.tFr_1BlkAllSk_Flag_Ct = tFr_Blk_Flag_AllSkOBwStdMean_Ct_Cell; 
oFlag_VctyEr_Out.OB.tFr_1BlkAllSk_Flag_Pct = tFr_Blk_Flag_AllSkOBwStdMean_Pct_Cell;

oFlag_VctyEr_Out.OB.tFrMidPt_1BlkAllSk_Flag_Ct = tFr_MidPt_Blk_Flag_AllSkOBwStdMean_Ct_Cell; 
oFlag_VctyEr_Out.OB.tFrMidPt_1BlkAllSk_Flag_Pct = tFr_MidPt_Blk_Flag_AllSkOBwStdMean_Pct_Cell;

oFlag_VctyEr_Out.OB.tMidPt_AllBlkAllSk_Flag_Ct = tMidPt_Flag_AllSkOBwStdMean_Ct; 
oFlag_VctyEr_Out.OB.tMidPt_AllBlkAllSk_Flag_Pct = tMidPt_Flag_AllSkOBwStdMean_Pct;

oFlag_VctyEr_Out.OB.tMidPt_AllBlkAllSk_Metric_Mean = tMidPt_Metric_AllSkOBwStdMean_Mean;
oFlag_VctyEr_Out.OB.tMidPt_AllBlkAllSk_Metric_Median = tMidPt_Metric_AllSkOBwStdMean_Median;

oFlag_VctyEr_Out.OB.tFr_AllBlkAllSk_Flag_Ct = tFr_Flag_AllSkOBwStdMean_Ct; 
oFlag_VctyEr_Out.OB.tFr_AllBlkAllSk_Flag_Pct = tFr_Flag_AllSkOBwStdMean_Pct;

oFlag_VctyEr_Out.OB.tFrMidPt_AllBlkAllSk_Flag_Ct = tFr_MidPt_Flag_AllSkOBwStdMean_Ct; 
oFlag_VctyEr_Out.OB.tFrMidPt_AllBlkAllSk_Flag_Pct = tFr_MidPt_Flag_AllSkOBwStdMean_Pct;

oFlag_VctyEr_Out.OB.tFrNonFlag_AllBlkAllSk_MaxContCt = MaxNonFlagtFrCt_AllSkOBwStdMean; 
oFlag_VctyEr_Out.OB.tFrNonFlag_AllBlkAllSk_MaxContCt_FrStart = MaxNonFlagtFrCt_AllSkOBwStdMean_FrStart;
oFlag_VctyEr_Out.OB.tFrNonFlag_AllBlkAllSk_MaxContCt_FrEnd = MaxNonFlagtFrCt_AllSkOBwStdMean_FrEnd;


%% Calculate Velocity error: ii. Velocity Error from Skeleton Statistics

tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell = cell2mat(oVelocity_In.tMidPt_1SkOBwMean_SkStdEr_Cell).*2./cell2mat(oVelocity_In.tMidPt_1SkOBwMean_SkMean_Cell); 
tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell = tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell > StdErLimSK_In; 

% If oVelocity.tMidPt_1SkOBct_Cell = 0 or 1; ie, if tMidPt doesn't
% have orient bands (and orient and velocity values are interpolated), also
% flag the cell. (include 1 as OB ct = 1yields errors in Stdev and StdEr
% values)
tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell = tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell | logical(max(cell2mat(oVelocity_In.tMidPt_1SkOBct_Cell)<2,[],1));

tMidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell = zeros(1,1,NumBlk_In,'single'); 
tMidPt_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell = zeros(1,1,NumBlk_In,'single'); 

tMidPtwMean_Blk_Metric_1SkOBwMean_SkMean_Cell = zeros(1,1,NumBlk_In,'single'); 
tMidPtwMedian_Blk_Metric_1SkOBwMean_SkMean_Cell = zeros(1,1,NumBlk_In,'single'); 

tFr_Blk_Flag_1SkOBwMean_SkMean_Cell = zeros(TimeParam_In.tMidPt_FrEnd(end),1,NumBlk_In,'logical'); 
tFr_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell = zeros(1,1,NumBlk_In,'single'); 
tFr_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell = zeros(1,1,NumBlk_In,'single'); 

for BlkCt = 1:NumBlk_In

    tMidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell(1,1,BlkCt) = sum(tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell(1,:,BlkCt),2);
    tMidPt_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell(1,1,BlkCt) = tMidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell(1,1,BlkCt)/NumMidPt_In*100;

    tMidPtwMean_Blk_Metric_1SkOBwMean_SkMean_Cell(1,1,BlkCt) = ...
        sum(TimeParam_In.tMidPt_FrLength.*tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell(1,:,BlkCt),2)/sum(TimeParam_In.tMidPt_FrLength,2);

    % weightedMedian() from FileExchange:
    % https://www.mathworks.com/matlabcentral/fileexchange/23077-weighted-median
    tMidPtwMedian_Blk_Metric_1SkOBwMean_SkMean_Cell(1,1,BlkCt) = ...
        weightedMedian(tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell(1,:,BlkCt),TimeParam_In.tMidPt_FrLength);

    tFr_Flag_Temp = zeros(TimeParam_In.tMidPt_FrEnd(end),1,1,'logical');
    tFr_Flag_Temp(ismember(tFr_HistBin,find(tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell(:,:,BlkCt))),1) = 1;
    tFr_Blk_Flag_1SkOBwMean_SkMean_Cell(:,1,BlkCt) = tFr_Flag_Temp;

    tFr_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell(1,1,BlkCt) = sum(tFr_HistCt(find(tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell(:,:,BlkCt))),2);
    tFr_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell(1,1,BlkCt) = tFr_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell(1,1,BlkCt)/TimeParam_In.tMidPt_FrEnd(end)*100;

end

tMidPt_Flag_1SkOBwMean_SkMean_Ct = sum(tMidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell,'all'); 
tMidPt_Flag_1SkOBwMean_SkMean_Pct = tMidPt_Flag_1SkOBwMean_SkMean_Ct/(NumBlk_In*NumMidPt_In)*100; 

tMidPt_Metric_1SkOBwMean_SkMean_Mean = mean(tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell,'all'); 
tMidPt_Metric_1SkOBwMean_SkMean_Median = median(tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell,'all'); 

tFr_Flag_1SkOBwMean_SkMean_Ct = sum(tFr_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell,'all'); 
tFr_Flag_1SkOBwMean_SkMean_Pct = tFr_Flag_1SkOBwMean_SkMean_Ct/(NumBlk_In*TimeParam_In.tMidPt_FrEnd(end))*100; 

% = Determine longest non-flagged stretch of time (for all vessel blocks)
rp_NonFlag_1SkOBwMean_SkMean = regionprops(~logical(max(tFr_Blk_Flag_1SkOBwMean_SkMean_Cell,[],3)),'Area','PixelIdxList'); 

if ~isempty(rp_NonFlag_1SkOBwMean_SkMean) % NonFlag regions exist
    [MaxNonFlagtFrCt_1SkOBwMean_SkMean, MaxNonFlagtFrCt_1SkOBwMean_SkMean_BlobIdx] = max(cat(1,rp_NonFlag_1SkOBwMean_SkMean.Area),[],'all');
    MaxNonFlagtFrCt_1SkOBwMean_SkMean_FrStart = rp_NonFlag_1SkOBwMean_SkMean(MaxNonFlagtFrCt_1SkOBwMean_SkMean_BlobIdx).PixelIdxList(1);
    MaxNonFlagtFrCt_1SkOBwMean_SkMean_FrEnd = rp_NonFlag_1SkOBwMean_SkMean(MaxNonFlagtFrCt_1SkOBwMean_SkMean_BlobIdx).PixelIdxList(end);
else 
    MaxNonFlagtFrCt_1SkOBwMean_SkMean = 0;
    MaxNonFlagtFrCt_1SkOBwMean_SkMean_FrStart = 0;
    MaxNonFlagtFrCt_1SkOBwMean_SkMean_FrEnd = 0;
end

% = Convert to cell
tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell = mat2cell(tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell,1,repelem(1,NumMidPt_In),repelem(1,NumBlk_In));
tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell = mat2cell(tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell,1,repelem(1,NumMidPt_In),repelem(1,NumBlk_In));
tMidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell = mat2cell(tMidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell,1,1,repelem(1,NumBlk_In));
tMidPt_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell = mat2cell(tMidPt_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell,1,1,repelem(1,NumBlk_In));

tMidPtwMean_Blk_Metric_1SkOBwMean_SkMean_Cell = mat2cell(tMidPtwMean_Blk_Metric_1SkOBwMean_SkMean_Cell,1,1,repelem(1,NumBlk_In));
tMidPtwMedian_Blk_Metric_1SkOBwMean_SkMean_Cell = mat2cell(tMidPtwMedian_Blk_Metric_1SkOBwMean_SkMean_Cell,1,1,repelem(1,NumBlk_In));

tFr_Blk_Flag_1SkOBwMean_SkMean_Cell = mat2cell(tFr_Blk_Flag_1SkOBwMean_SkMean_Cell,TimeParam_In.tMidPt_FrEnd(end),1,repelem(1,NumBlk_In));
tFr_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell = mat2cell(tFr_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell,1,1,repelem(1,NumBlk_In));
tFr_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell = mat2cell(tFr_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell,1,1,repelem(1,NumBlk_In));

tFr_MidPt_Blk_Flag_1SkOBwMean_SkMean_Cell = tFr_Blk_Flag_1SkOBwMean_SkMean_Cell;
tFr_MidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell = tFr_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell;
tFr_MidPt_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell = tFr_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell;
tFr_MidPt_Flag_1SkOBwMean_SkMean_Ct = tFr_Flag_1SkOBwMean_SkMean_Ct;
tFr_MidPt_Flag_1SkOBwMean_SkMean_Pct = tFr_Flag_1SkOBwMean_SkMean_Pct; 

% = Load into struct
oFlag_VctyEr_Out.SK.tMidPt_1BlkAllSk_Flag_Cell = tMidPt_Blk_Flag_1SkOBwMean_SkMean_Cell; 
oFlag_VctyEr_Out.SK.tMidPt_1BlkAllSk_Metric_Cell = tMidPt_Blk_Metric_1SkOBwMean_SkMean_Cell; 

oFlag_VctyEr_Out.SK.tFr_1BlkAllSk_Flag_Cell = tFr_Blk_Flag_1SkOBwMean_SkMean_Cell; 

oFlag_VctyEr_Out.SK.tMidPt_1BlkAllSk_Flag_Ct = tMidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell;  
oFlag_VctyEr_Out.SK.tMidPt_1BlkAllSk_Flag_Pct = tMidPt_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell;

oFlag_VctyEr_Out.SK.tMidPtwMean_1BlkAllSk_Metric_Cell = tMidPtwMean_Blk_Metric_1SkOBwMean_SkMean_Cell;
oFlag_VctyEr_Out.SK.tMidPtwMedian_1BlkAllSk_Metric_Cell = tMidPtwMedian_Blk_Metric_1SkOBwMean_SkMean_Cell;

oFlag_VctyEr_Out.SK.tFr_1BlkAllSk_Flag_Ct = tFr_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell; 
oFlag_VctyEr_Out.SK.tFr_1BlkAllSk_Flag_Pct = tFr_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell;

oFlag_VctyEr_Out.SK.tFrMidPt_1BlkAllSk_Flag_Ct = tFr_MidPt_Blk_Flag_1SkOBwMean_SkMean_Ct_Cell; 
oFlag_VctyEr_Out.SK.tFrMidPt_1BlkAllSk_Flag_Pct = tFr_MidPt_Blk_Flag_1SkOBwMean_SkMean_Pct_Cell;

oFlag_VctyEr_Out.SK.tMidPt_AllBlkAllSk_Flag_Ct = tMidPt_Flag_1SkOBwMean_SkMean_Ct; 
oFlag_VctyEr_Out.SK.tMidPt_AllBlkAllSk_Flag_Pct = tMidPt_Flag_1SkOBwMean_SkMean_Pct;

oFlag_VctyEr_Out.SK.tMidPt_AllBlkAllSk_Metric_Mean = tMidPt_Metric_1SkOBwMean_SkMean_Mean;
oFlag_VctyEr_Out.SK.tMidPt_AllBlkAllSk_Metric_Median = tMidPt_Metric_1SkOBwMean_SkMean_Median;

oFlag_VctyEr_Out.SK.tFr_AllBlkAllSk_Flag_Ct = tFr_Flag_1SkOBwMean_SkMean_Ct; 
oFlag_VctyEr_Out.SK.tFr_AllBlkAllSk_Flag_Pct = tFr_Flag_1SkOBwMean_SkMean_Pct;

oFlag_VctyEr_Out.SK.tFrMidPt_AllBlkAllSk_Flag_Ct = tFr_MidPt_Flag_1SkOBwMean_SkMean_Ct; 
oFlag_VctyEr_Out.SK.tFrMidPt_AllBlkAllSk_Flag_Pct = tFr_MidPt_Flag_1SkOBwMean_SkMean_Pct;

oFlag_VctyEr_Out.SK.tFrNonFlag_AllBlkAllSk_MaxContCt = MaxNonFlagtFrCt_1SkOBwMean_SkMean; 
oFlag_VctyEr_Out.SK.tFrNonFlag_AllBlkAllSk_MaxContCt_FrStart = MaxNonFlagtFrCt_1SkOBwMean_SkMean_FrStart;
oFlag_VctyEr_Out.SK.tFrNonFlag_AllBlkAllSk_MaxContCt_FrEnd = MaxNonFlagtFrCt_1SkOBwMean_SkMean_FrEnd;


%% Report on screen

% % % fprintf(strcat('\n====================\n')); % Must match fCleanImgFFTFilteredCell() value
% % % 
% % % for vmCt = 1:2
% % % 
% % %     if vmCt == 1 % variation method
% % %         s = oFlag_VctyEr_Out.OB;
% % %     elseif vmCt == 2
% % %         s = oFlag_VctyEr_Out.SK;
% % %     end
% % % 
% % %     fprintf(strcat(num2str(s.tMidPt_1BlkAllSk_Flag_Ct{1,1,ScrnReportBlkChoice_In},'%d'),'/',num2str(NumMidPt_In,'%d'),32,'(',num2str(s.tMidPt_1BlkAllSk_Flag_Pct{1,1,ScrnReportBlkChoice_In},'%0.2f'),'%%',')',32,'tMidPt flagged for large velocity standard error.\n')); 
% % %     fprintf(strcat(num2str(s.tFr_1BlkAllSk_Flag_Ct{1,1,ScrnReportBlkChoice_In},'%d'),'/',num2str(TimeParam_In.tSeg_FrEnd(end),'%d'),32,'(',num2str(s.tFr_1BlkAllSk_Flag_Pct{1,1,ScrnReportBlkChoice_In},'%0.2f'),'%%',')',32,'tFr flagged for large velocity standard error.\n')); 
% % %     fprintf(strcat(num2str(s.tFrMidPt_1BlkAllSk_Flag_Ct{1,1,ScrnReportBlkChoice_In},'%d'),'/',num2str(TimeParam_In.tSeg_FrEnd(end),'%d'),32,'(',num2str(s.tFrMidPt_1BlkAllSk_Flag_Pct{1,1,ScrnReportBlkChoice_In},'%0.2f'),'%%',')',32,'tFr in tMidPts flagged for large velocity standard error.\n'));
% % % 
% % %     % Issue warning if either % flagged tMidPt and % flagged tFr is greater than 20%
% % %     if s.tMidPt_1BlkAllSk_Flag_Pct{1,1,ScrnReportBlkChoice_In} > 20 | s.tFr_1BlkAllSk_Flag_Pct{1,1,ScrnReportBlkChoice_In} > 20
% % %         warning('***Significant unstable flow may be present in input video. Velocity measurements may not be representative of true values.***');
% % %     end
% % % 
% % %     fprintf('\n');
% % % 
% % % fprintf(strcat('Longest stretch of low velocity standard error time frames (all vessel blocks considered):\n'));
% % % fprintf(strcat('Time frame count:',32,num2str(s.tFrNonFlag_AllBlkAllSk_MaxContCt,'%d'),'\n'));
% % % fprintf(strcat('Start frame:',32,num2str(s.tFrNonFlag_AllBlkAllSk_MaxContCt_FrStart,'%d'),';',32,'End frame:',32,num2str(s.tFrNonFlag_AllBlkAllSk_MaxContCt_FrEnd,'%d'),'\n'));
% % % 
% % % end
% % % 
% % % fprintf(strcat('====================\n')); 


