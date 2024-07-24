function [oFlag_Out] = fConvertFlag_tFrtoMidPtAllSk(tFr_Blk_Flag_In,Frac_tMidPtFrLength_Flag_In, TimeParam_In)

NumSeg_In = TimeParam_In.Num_tSeg;   
NumMidPt_In = TimeParam_In.Num_tMidPt;
NumBlk_In = size(tFr_Blk_Flag_In,3); 

tFr_HistEdge = horzcat(TimeParam_In.tMidPt_FrStart,TimeParam_In.tMidPt_FrEnd(end));
[tFr_HistCt,~,tFr_HistBin] = histcounts(1:1:TimeParam_In.tMidPt_FrEnd(end),tFr_HistEdge);


%% Create Flag_Cell of {SkelCt,MidPtCt,BlkCt} format based on fraction of time frames likely w/ Registration Errors

if islogical(tFr_Blk_Flag_In) & tFr_Blk_Flag_In(1)>-9999  

    fprintf('Creating Flag_Cell based on fraction of flagged tFr...\n');
    
    tFr_Blk_Flag_Cell = mat2cell(tFr_Blk_Flag_In,TimeParam_In.tMidPt_FrEnd(end),1,repelem(1,NumBlk_In)); 

    tMidPt_Blk_Flag_Cell = cell(1,NumMidPt_In,NumBlk_In); 
    tFr_MidPt_Blk_Flag_Cell = repmat({zeros(TimeParam_In.tSeg_FrEnd(end),1,'logical')},1,1,NumBlk_In); 

    tMidPt_Blk_Flag_Ct_Cell = cell(1,1,NumBlk_In); 
    tMidPt_Blk_Flag_Pct_Cell = cell(1,1,NumBlk_In); 

    tFr_Blk_Flag_Ct_Cell = cell(1,1,NumBlk_In); 
    tFr_Blk_Flag_Pct_Cell = cell(1,1,NumBlk_In); 

    tFr_MidPt_Blk_Flag_Ct_Cell = cell(1,1,NumBlk_In); 
    tFr_MidPt_Blk_Flag_Pct_Cell = cell(1,1,NumBlk_In); 

    for BlkCt = 1:NumBlk_In

        rp = regionprops(tFr_Blk_Flag_In(:,:,BlkCt),'PixelIdxList');
        tFr_Blk_List = cat(1,rp.PixelIdxList); 

        for MidPtCt = 1:NumMidPt_In

            tFr_Blk_inMidPt_Idx = ismember(tFr_Blk_List,TimeParam_In.tMidPt_FrStart(MidPtCt):TimeParam_In.tMidPt_FrEnd(MidPtCt)); 

            if sum(tFr_Blk_inMidPt_Idx,'all')/(tFr_HistCt(MidPtCt))>Frac_tMidPtFrLength_Flag_In
                tMidPt_Blk_Flag_Cell{1,MidPtCt,BlkCt} = logical(1);
                tFr_MidPt_Blk_Flag_Cell{1,1,BlkCt}(TimeParam_In.tMidPt_FrStart(MidPtCt):TimeParam_In.tMidPt_FrEnd(MidPtCt),1) = 1;
            else
                tMidPt_Blk_Flag_Cell{1,MidPtCt,BlkCt} = logical(0);
            end

        end

        tMidPt_Blk_Flag_Ct_Cell{1,1,BlkCt} = sum(cell2mat(tMidPt_Blk_Flag_Cell(1,:,BlkCt)),'all');
        tMidPt_Blk_Flag_Pct_Cell{1,1,BlkCt} = tMidPt_Blk_Flag_Ct_Cell{1,1,BlkCt}/NumMidPt_In*100;

        tFr_Blk_Flag_Ct_Cell{1,1,BlkCt} = sum(cell2mat(tFr_Blk_Flag_Cell(1,:,BlkCt)),'all');
        tFr_Blk_Flag_Pct_Cell{1,1,BlkCt} = tFr_Blk_Flag_Ct_Cell{1,1,BlkCt}/TimeParam_In.tMidPt_FrEnd(end)*100;

        tFr_MidPt_Blk_Flag_Ct_Cell{1,1,BlkCt} = sum(cell2mat(tFr_MidPt_Blk_Flag_Cell(1,:,BlkCt)),'all');
        tFr_MidPt_Blk_Flag_Pct_Cell{1,1,BlkCt} = tFr_MidPt_Blk_Flag_Ct_Cell{1,1,BlkCt}/TimeParam_In.tMidPt_FrEnd(end)*100;

    end

    tMidPt_Flag_Ct = sum(cell2mat(tMidPt_Blk_Flag_Ct_Cell),'all'); 
    tMidPt_Flag_Pct = tMidPt_Flag_Ct/(NumBlk_In*NumMidPt_In)*100;

    tFr_Flag_Ct = sum(cell2mat(tFr_Blk_Flag_Ct_Cell),'all'); 
    tFr_Flag_Pct = tFr_Flag_Ct/(NumBlk_In*TimeParam_In.tMidPt_FrEnd(end))*100;

    tFr_MidPt_Flag_Ct = sum(cell2mat(tFr_MidPt_Blk_Flag_Ct_Cell),'all'); 
    tFr_MidPt_Flag_Pct = tFr_MidPt_Flag_Ct/(NumBlk_In*TimeParam_In.tMidPt_FrEnd(end))*100;

else

    tMidPt_Blk_Flag_Cell = repmat({-9999},1,NumMidPt_In,NumBlk_In);
    tFr_Blk_Flag_Cell = repmat({tFr_Blk_Flag_In},1,1,NumBlk_In);
    tFr_MidPt_Blk_Flag_Cell = repmat({tFr_Blk_Flag_In},1,1,NumBlk_In);

    tMidPt_Blk_Flag_Ct_Cell = repmat({-9999},1,1,NumBlk_In);
    tMidPt_Blk_Flag_Pct_Cell = repmat({-9999},1,1,NumBlk_In);
    tFr_Blk_Flag_Ct_Cell = repmat({-9999},1,1,NumBlk_In);
    tFr_Blk_Flag_Pct_Cell = repmat({-9999},1,1,NumBlk_In);
    tFr_MidPt_Blk_Flag_Ct_Cell = repmat({-9999},1,1,NumBlk_In);
    tFr_MidPt_Blk_Flag_Ct_Cell = repmat({-9999},1,1,NumBlk_In);

    tMidPt_Flag_Ct = -9999;
    tMidPt_Flag_Pct = -9999; 
    tFr_Flag_Ct = -9999;
    tFr_Flag_Pct = -9999;
    tFr_MidPt_Flag_Ct = -9999;
    tFr_MidPt_Flag_Pct = -9999;

end


%% Locate longest stretch of tFr without flag

if islogical(tFr_Blk_Flag_In) & tFr_Blk_Flag_In(1)>-9999 

    rp_NonFlag = regionprops(~logical(max(cell2mat(tFr_Blk_Flag_Cell),[],3)),'Area','PixelIdxList');

    if ~isempty(rp_NonFlag) 
        [MaxNonFlagtFrCt, MaxNonFlagtFrCt_BlobIdx] = max(cat(1,rp_NonFlag.Area),[],'all');
        MaxNonFlagtFrCt_FrStart = rp_NonFlag(MaxNonFlagtFrCt_BlobIdx).PixelIdxList(1);
        MaxNonFlagtFrCt_FrEnd = rp_NonFlag(MaxNonFlagtFrCt_BlobIdx).PixelIdxList(end);
    else 
        MaxNonFlagtFrCt = 0;
        MaxNonFlagtFrCt_FrStart = 0;
        MaxNonFlagtFrCt_FrEnd = 0;
    end

else

    MaxNonFlagtFrCt = -9999;
    MaxNonFlagtFrCt_FrStart = -9999;
    MaxNonFlagtFrCt_FrEnd = -9999;

end

clearvars rp_NonFlag MaxNonFlagtFrCt_BlobIdx;


%% Prepare Output Structure

oFlag_Out = struct;

oFlag_Out.tMidPt_1BlkAllSk_Flag_Cell = tMidPt_Blk_Flag_Cell; 
oFlag_Out.tFr_1BlkAllSk_Flag_Cell = tFr_Blk_Flag_Cell; 

oFlag_Out.tMidPt_1BlkAllSk_Flag_Ct = tMidPt_Blk_Flag_Ct_Cell; 
oFlag_Out.tMidPt_1BlkAllSk_Flag_Pct = tMidPt_Blk_Flag_Pct_Cell;

oFlag_Out.tFr_1BlkAllSk_Flag_Ct = tFr_Blk_Flag_Ct_Cell; 
oFlag_Out.tFr_1BlkAllSk_Flag_Pct = tFr_Blk_Flag_Pct_Cell;

oFlag_Out.tFrMidPt_1BlkAllSk_Flag_Ct = tFr_MidPt_Blk_Flag_Ct_Cell; 
oFlag_Out.tFrMidPt_1BlkAllSk_Flag_Pct = tFr_MidPt_Blk_Flag_Pct_Cell;

oFlag_Out.tMidPt_AllBlkAllSk_Flag_Ct = tMidPt_Flag_Ct; 
oFlag_Out.tMidPt_AllBlkAllSk_Flag_Pct = tMidPt_Flag_Pct; 

oFlag_Out.tFr_AllBlkAllSk_Flag_Ct = tFr_Flag_Ct;
oFlag_Out.tFr_AllBlkAllSk_Flag_Pct = tFr_Flag_Pct;

oFlag_Out.tFrMidPt_AllBlkAllSk_FlagCt = tFr_MidPt_Flag_Ct; 
oFlag_Out.tFrMidPt_AllBlkAllSk_FlagPct = tFr_MidPt_Flag_Pct;

oFlag_Out.tFrNonFlag_AllBlkAllSk_MaxContCt = MaxNonFlagtFrCt; 
oFlag_Out.tFrNonFlag_AllBlkAllSk_MaxContCt_FrStart = MaxNonFlagtFrCt_FrStart;
oFlag_Out.tFrNonFlag_AllBlkAllSk_MaxContCt_FrEnd = MaxNonFlagtFrCt_FrEnd;



