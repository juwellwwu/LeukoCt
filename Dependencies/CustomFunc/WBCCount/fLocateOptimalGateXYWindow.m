function [GateXYParam_Out,Fig_SkelPxStartIdx_BestWin_Out] = ...
    fLocateOptimalGateXYWindow(ImgStack_In,Img_FocusDiameter_Idx_In,...
    WinSkelPxIdx_Cell_In,WBCtoGateXYWinAreaRatio_In,...
    oSCS_In,...
    oRadius_In,...
    XY_PxLength_In,...
    GateXY_ByBlk_Option_In,... 
    GateXY_SingleWinLength_Option_In,... 
    GateXY_SingleWinLengthUserChoice_Option_In, GateXY_SingleWinLengthUserChoice_Idx_In) 


Skel_xy_Cell_In = oSCS_In.Skel_xy_Cell;
Skel_bdDist_Cell_In = oSCS_In.Skel_bdDist_Cell;

CellDiameter_px_In = 5/XY_PxLength_In;
NumBlk_In = size(Skel_xy_Cell_In,3);
SkelBlockLength_In = height(Skel_xy_Cell_In{1,1,1});

WBCArea_px = pi*(0.5*12.0/XY_PxLength_In).^2;

tMidPt_fDiameter_In = oRadius_In.F.tMidPt_SkPx_Diameter; 
tMidPt_lDiameter_In = oRadius_In.L.tMidPt_SkPx_Diameter; 
tMidPt_lDiameter_AllSk_Cell_In = oRadius_In.L.tMidPt_AllSk_Diameter_Cell; 


%% Re-organize Skel (x,y) 

[Skel_LinIdx_Mtx,~,~] = fLinearizeSkelxyCell(ImgStack_In,Skel_xy_Cell_In);


%% Calculate optimal window

if GateXY_ByBlk_Option_In == 1 % [Removed]

elseif GateXY_ByBlk_Option_In == 0 

    %% METHOD 4: Determine to be Flagged SkelPxCt with very low Img_FocusDIameter_Idx

    NumSkel_Temp = floor(0.5*0.8*width(Skel_LinIdx_Mtx))*2+1; 

    SkelPx_FocusDiameterIdxMean = zeros(height(Skel_LinIdx_Mtx),1,'single');
    for SkelPxCt = 1:height(Skel_LinIdx_Mtx)
        SkelPx_FocusDiameterIdxMean(SkelPxCt,1) = ...
            mean(Img_FocusDiameter_Idx_In(Skel_LinIdx_Mtx(SkelPxCt,1:NumSkel_Temp)),2);
    end
    
    % By stretch of SkelPx w/ consistently low focus
    % https://www.mathworks.com/matlabcentral/answers/253093-how-to-find-all-blocks-of-numbers-with-consecutive-locations-which-are-less-than-a-value-in-an-array
    threshold = 0.50*max(SkelPx_FocusDiameterIdxMean,[],'all'); % Default: 0.50
    transitions = diff([0, SkelPx_FocusDiameterIdxMean' < threshold, 0]);
    runstarts = find(transitions == 1);
    runends = find(transitions == -1) - 1;
    blocks = arrayfun(@(s, e) SkelPx_FocusDiameterIdxMean(s:e), runstarts, runends, 'UniformOutput', false);
    blockheight = cellfun(@height, blocks);
    
    height_threshold = floor(CellDiameter_px_In);
    SkelPx_FocusDiameterIdxMean_FlagID = zeros(size(SkelPx_FocusDiameterIdxMean),'logical');
    runstarts = runstarts(blockheight>height_threshold); 
    runends = runends(blockheight>height_threshold); 
    for runCt = 1:numel(runstarts)
        SkelPx_FocusDiameterIdxMean_FlagID(runstarts(runCt):runends(runCt)) = 1;
    end
    SkelPx_FocusDiameterIdxMean_FlagID = find(SkelPx_FocusDiameterIdxMean_FlagID);


    %% METHOD 4: Designate Window for Best Focus Location

     % = Number of WinArea to try       
     Num_WinArea_In = width(WinSkelPxIdx_Cell_In); % # Win Area tested (per SkelPxCt)
     
     % = Prepare cells that hold all window information for each WinArea
     Num_Win =  height(Skel_LinIdx_Mtx); 

     FocusDiameterIdxMean_AllWin_Cell = repmat({zeros(Num_Win,1,'single')},1,Num_WinArea_In); 
     VesselDiameter_px_AllWin_Cell = repmat({zeros(Num_Win,1,'single')},1,Num_WinArea_In); 
     NumSkel_Allwin_Cell = repmat({zeros(Num_Win,1,'single')},1,Num_WinArea_In); 

     FocusDiameterIdxMean_BestWin_SkelPxStartIdx = zeros(1,Num_WinArea_In,'single'); 
     NumSkel_BestWin = zeros(1,Num_WinArea_In,'single'); 
     VesselDiameter_px_BestWin = zeros(1,Num_WinArea_In,'single'); 
     WinLength_px_BestWin = zeros(1,Num_WinArea_In,'single'); 
     WinSkelPxLinIdx_BestWin_Cell = cell(1,Num_WinArea_In);

     FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell = cell(1,Num_WinArea_In); 
     NumSkel_BestMultiWin_Cell = cell(1,Num_WinArea_In);
     VesselDiameter_px_BestMultiWin_Cell = cell(1,Num_WinArea_In);
     WinLength_px_BestMultiWin_Cell = cell(1,Num_WinArea_In);
     WinSkelPxLinIdx_BestMultiWin_Cell = cell(1,Num_WinArea_In);

     for aCt = 1:Num_WinArea_In 

        for wCt = 1:Num_Win 

            NumSkel_Allwin_Cell{1,aCt}(wCt,1) = (floor(0.5*0.8*width(Skel_LinIdx_Mtx))*2+1);

            % = Which skeleton pixels to include in window
            if ~isnan(WinSkelPxIdx_Cell_In{1,aCt}(wCt,2))

                Win_SkelPxStartIdx = WinSkelPxIdx_Cell_In{1,aCt}(wCt,1);
                Win_SkelPxEndIdx = WinSkelPxIdx_Cell_In{1,aCt}(wCt,2);

                Win_LinIdx = Skel_LinIdx_Mtx(Win_SkelPxStartIdx:Win_SkelPxEndIdx,1:NumSkel_Allwin_Cell{1,aCt}(wCt,1));
                Win_LinIdx = reshape(Win_LinIdx,numel(Win_LinIdx),1);

                % = Quality of focus + diameter in window (better = higher value)
                if any(ismember(Win_SkelPxStartIdx:Win_SkelPxEndIdx,SkelPx_FocusDiameterIdxMean_FlagID))
                    FocusDiameterIdxMean_AllWin_Cell{1,aCt}(wCt,1) = NaN;
                else
                    FocusDiameterIdxMean_AllWin_Cell{1,aCt}(wCt,1) = mean(Img_FocusDiameter_Idx_In(Win_LinIdx),'all');
                end

                % = Mean vessel diameter within window
                VesselDiameter_px_AllWin_Cell{1,aCt}(wCt,1) = mean(max(tMidPt_lDiameter_In(Win_SkelPxStartIdx:Win_SkelPxEndIdx,:),[],2),1);
                
            else

                FocusDiameterIdxMean_AllWin_Cell{1,aCt}(wCt,1) = NaN;
                VesselDiameter_px_AllWin_Cell{1,aCt}(wCt,1) = NaN;

            end

        end

        % = Find window index of Best focus + diameter Window for each WinArea
        [FocusDiameterIdxMean_AllWin_Sort,FocusDiameterIdxMean_AllWin_SortIdx] = sort(FocusDiameterIdxMean_AllWin_Cell{1,aCt},1,'descend','MissingPlacement','last');

        FocusDiameterIdxMean_AllWin_SortIdx(isnan(FocusDiameterIdxMean_AllWin_Sort)) = []; 
        FocusDiameterIdxMean_AllWin_Sort(isnan(FocusDiameterIdxMean_AllWin_Sort)) = [];

        FocusDiameterIdxMean_AllWin_SortIdx(FocusDiameterIdxMean_AllWin_Sort<0.75*max(FocusDiameterIdxMean_AllWin_Sort,[],'all')) = []; 
        FocusDiameterIdxMean_AllWin_Sort(FocusDiameterIdxMean_AllWin_Sort<0.75*max(FocusDiameterIdxMean_AllWin_Sort,[],'all')) = [];

        if ~isempty(FocusDiameterIdxMean_AllWin_SortIdx)

            % i) Identify best "1" window index that qualify for cell counting
            FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt) = FocusDiameterIdxMean_AllWin_SortIdx(1,1); 
            NumSkel_BestWin(1,aCt) = NumSkel_Allwin_Cell{1,aCt}(FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt),1); 
            VesselDiameter_px_BestWin(1,aCt) = VesselDiameter_px_AllWin_Cell{1,aCt}(FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt),1); 
            WinLength_px_BestWin(1,aCt) = ... 
                WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt),2)-WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt),1)+1;
            WinSkelPxLinIdx_BestWin_Cell{1,aCt} = ...
                Skel_LinIdx_Mtx(WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt),1):WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt),2),1:NumSkel_BestWin(1,aCt));
            WinSkelPxLinIdx_BestWin_Cell{1,aCt} = reshape(WinSkelPxLinIdx_BestWin_Cell{1,aCt}, numel(WinSkelPxLinIdx_BestWin_Cell{1,aCt}),1);

            % ii) Identify "all" window indices that qualify for cell counting
            FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt} = FocusDiameterIdxMean_AllWin_SortIdx(1,1); 
            WinSkelPx_Temp = ... 
                WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_AllWin_SortIdx(1,1),1): WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_AllWin_SortIdx(1,1),2); 
            
            for iCt = 2:numel(FocusDiameterIdxMean_AllWin_SortIdx) 
                SkelPx_StartIdx_Temp = WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_AllWin_SortIdx(iCt,1),1);
                SkelPx_EndIdx_Temp = WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_AllWin_SortIdx(iCt,1),2);
                if ~ismember(SkelPx_StartIdx_Temp:SkelPx_EndIdx_Temp,WinSkelPx_Temp)
                    FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt} = vertcat(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt},FocusDiameterIdxMean_AllWin_SortIdx(iCt,1));
                    WinSkelPx_Temp = union(WinSkelPx_Temp,SkelPx_StartIdx_Temp:SkelPx_EndIdx_Temp); 
                end
            end

            NumSkel_BestMultiWin_Cell{1,aCt} = NumSkel_Allwin_Cell{1,aCt}(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt},1); 
            VesselDiameter_px_BestMultiWin_Cell{1,aCt} = VesselDiameter_px_AllWin_Cell{1,aCt}(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt},1); 
            WinLength_px_BestMultiWin_Cell{1,aCt} = ...
                WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt},2)-WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt},1)+1;
            
            WinSkelPxLinIdx_BestMultiWin_Cell{1,aCt} = cell(numel(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt}),1); 
            for iCt = 1:numel(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt})
                WinSkelPxLinIdx_BestMultiWin_Cell{1,aCt}{iCt,1} = ...
                    Skel_LinIdx_Mtx(WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt}(iCt,1),1):WinSkelPxIdx_Cell_In{1,aCt}(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt}(iCt,1),2),1:NumSkel_BestMultiWin_Cell{1,aCt}(iCt,1));
                WinSkelPxLinIdx_BestMultiWin_Cell{1,aCt}{iCt,1} = reshape(WinSkelPxLinIdx_BestMultiWin_Cell{1,aCt}{iCt,1}, numel(WinSkelPxLinIdx_BestMultiWin_Cell{1,aCt}{iCt,1}),1);
            end

        else

            FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt) = NaN;
            NumSkel_BestWin(1,aCt) = NaN;
            VesselDiameter_px_BestWin(1,aCt) = NaN;
            WinLength_px_BestWin(1,aCt) = NaN;
            WinSkelPxLinIdx_BestWin_Cell{1,aCt} = NaN;

            FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,aCt} = NaN;
            NumSkel_BestMultiWin_Cell{1,aCt} = NaN;
            VesselDiameter_px_BestMultiWin_Cell{1,aCt} = NaN; 
            WinLength_px_BestMultiWin_Cell{1,aCt} = NaN;
            WinSkelPxLinIdx_BestMultiWin_Cell{1,aCt} = NaN;

        end

    end
    
    % = Plot: Multiple WinArea BestWin SkelPxStartIdx
    FigHandle = figure('Position',[150 600 450 500],'visible','off');
    cmap = colormap(copper(round(Num_WinArea_In*1))); 
    for aCt = 1:Num_WinArea_In 
        hold on
        plot(FocusDiameterIdxMean_AllWin_Cell{1,aCt},'Color',cmap(aCt,:),'Linewidth',1);
        xline(FocusDiameterIdxMean_BestWin_SkelPxStartIdx(1,aCt),'--','Color',cmap(aCt,:));
    end
    xlabel('Win_SkelPxStartIdx','Interpreter','none');
    ylabel('Win_FocusDiameterIdxMean','Interpreter','none');
    title({strcat('Win_FocusDiameterIdxMean: Best Window SkelPxStartIdx =',32,strjoin(arrayfun(@num2str,FocusDiameterIdxMean_BestWin_SkelPxStartIdx,'UniformOutput', 0),', ')),...
        strcat('Mean value, Win_FocusDiameterIdxMean =',32,strjoin(cellfun(@(x) num2str(mean(x,'all','omitmissing'),'%0.2f'),FocusDiameterIdxMean_AllWin_Cell,'UniformOutput',false),', ')),...
        strcat('Stdev value, Win_FocusDiameterIdxMean =',32,strjoin(cellfun(@(x) num2str(std(x,0,'all','omitmissing'),'%0.2f'),FocusDiameterIdxMean_AllWin_Cell,'UniformOutput',false),', ')),...
        strcat('WBCtoGateXYWinAreaRatio =',32,strjoin(cellfun(@(x) num2str(x),mat2cell(WBCtoGateXYWinAreaRatio_In,1,repelem(1,numel(WBCtoGateXYWinAreaRatio_In))),'UniformOutput',false),', ')),...
        strcat('WinLength (um) =',32,strjoin(cellfun(@(x) num2str(x,'%0.2f'),mat2cell(WinLength_px_BestWin.*XY_PxLength_In,1,repelem(1,numel(WinLength_px_BestWin))),'UniformOutput',false),', ')),...
        strcat('# NaN (low focus Win_StartPxIdx) =',32,strjoin(cellfun(@(x) num2str(numel(find(isnan(x))),'%d'),FocusDiameterIdxMean_AllWin_Cell,'UniformOutput',false),', '),32,'(',strjoin(cellfun(@(x) num2str(numel(find(isnan(x)))/height(x)*100,'%0.2f'),FocusDiameterIdxMean_AllWin_Cell,'UniformOutput',false),', '),'%)')},...
        'Interpreter','none');
    fontsize(7,'points');

    FigHandle_fr=getframe(gcf);
    Fig_SkelPxStartIdx_BestWin_Out = FigHandle_fr.cdata;

end


%% Prepare Output 

% = HistEdge for sorting SkelPx into VesselBlk bins
SkelPx_Blk_HistEdge = horzcat(1:SkelBlockLength_In:height(Skel_LinIdx_Mtx),height(Skel_LinIdx_Mtx));

if GateXY_ByBlk_Option_In == 1 % [Removed]

elseif GateXY_ByBlk_Option_In == 0 

    if GateXY_SingleWinLength_Option_In == 0

        if any(~isnan(FocusDiameterIdxMean_BestWin_SkelPxStartIdx))

            GateXYParam_Out.WinLength_px = WinLength_px_BestWin;
            GateXYParam_Out.WBCtoGateXYWinAreaRatio = WBCtoGateXYWinAreaRatio_In;
            GateXYParam_Out.WinVesselDiameter_px = VesselDiameter_px_BestWin;
            GateXYParam_Out.WinArea_px = GateXYParam_Out.WinVesselDiameter_px.*GateXYParam_Out.WinLength_px; 
            GateXYParam_Out.WinSkelPxStartIdx = FocusDiameterIdxMean_BestWin_SkelPxStartIdx;
            GateXYParam_Out.WinSkelPxBlk = nan(1,numel(GateXYParam_Out.WinSkelPxStartIdx),'single'); 
            for BlkCt = 1:numel(GateXYParam_Out.WinSkelPxStartIdx)
                if ~isnan(GateXYParam_Out.WinSkelPxStartIdx(1,BlkCt)) & ~isnan(GateXYParam_Out.WinLength_px(1,BlkCt))
                    [SkelPx_Blk_HistCt] = histcounts((GateXYParam_Out.WinSkelPxStartIdx(1,BlkCt):GateXYParam_Out.WinSkelPxStartIdx(1,BlkCt)+GateXYParam_Out.WinLength_px(1,BlkCt)),...
                        SkelPx_Blk_HistEdge);
                    [~,GateXYParam_Out.WinSkelPxBlk(1,BlkCt)] = max(SkelPx_Blk_HistCt);
                end
            end
            GateXYParam_Out.WinSkelPxLinIdx_Cell = WinSkelPxLinIdx_BestWin_Cell; 
            
            GateXYParam_Out.WinLength_px(isnan(WinLength_px_BestWin)) = [];
            GateXYParam_Out.WBCtoGateXYWinAreaRatio(isnan(WinLength_px_BestWin)) = [];
            GateXYParam_Out.WinVesselDiameter_px(isnan(WinLength_px_BestWin)) = [];
            GateXYParam_Out.WinArea_px(isnan(WinLength_px_BestWin)) = [];
            GateXYParam_Out.WinSkelPxStartIdx(isnan(WinLength_px_BestWin)) = [];
            GateXYParam_Out.WinSkelPxBlk(isnan(WinLength_px_BestWin)) = [];
            GateXYParam_Out.WinSkelPxLinIdx_Cell(isnan(WinLength_px_BestWin)) = [];

        else % if all window areas have no best window

            GateXYParam_Out.WinSkelPxStartIdx = NaN;
            GateXYParam_Out.WinVesselDiameter_px = NaN;
            warning('No Window SkelPxStartIdx found that meets focus requirements.')

        end

    elseif GateXY_SingleWinLength_Option_In == 1

        if GateXY_SingleWinLengthUserChoice_Option_In == 0 
            [~,SingleWinLengthIdx] = ... 
                max(cellfun(@height,FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell),[],2);
        elseif GateXY_SingleWinLengthUserChoice_Option_In == 1
            SingleWinLengthIdx = GateXY_SingleWinLengthUserChoice_Idx_In;
        end

        if ~isnan(FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,SingleWinLengthIdx}(1,1))
            GateXYParam_Out.WinLength_px = (WinLength_px_BestMultiWin_Cell{1,SingleWinLengthIdx})';
            GateXYParam_Out.WBCtoGateXYWinAreaRatio = repmat(WBCtoGateXYWinAreaRatio_In(1,SingleWinLengthIdx),1,numel(GateXYParam_Out.WinLength_px));
            GateXYParam_Out.WinVesselDiameter_px = (VesselDiameter_px_BestMultiWin_Cell{1,SingleWinLengthIdx})';
            GateXYParam_Out.WinArea_px =  GateXYParam_Out.WinVesselDiameter_px.*GateXYParam_Out.WinLength_px; 
            GateXYParam_Out.WinSkelPxStartIdx = (FocusDiameterIdxMean_BestMultiWin_SkelPxStartIdx_Cell{1,SingleWinLengthIdx})';
            GateXYParam_Out.WinSkelPxBlk = nan(1,numel(GateXYParam_Out.WinSkelPxStartIdx),'single');
            for BlkCt = 1:numel(GateXYParam_Out.WinSkelPxStartIdx)
                if ~isnan(GateXYParam_Out.WinSkelPxStartIdx(1,BlkCt)) & ~isnan(GateXYParam_Out.WinLength_px(1,BlkCt))
                    [SkelPx_Blk_HistCt] = histcounts((GateXYParam_Out.WinSkelPxStartIdx(1,BlkCt):GateXYParam_Out.WinSkelPxStartIdx(1,BlkCt)+GateXYParam_Out.WinLength_px(1,BlkCt)),...
                        SkelPx_Blk_HistEdge);
                    [~,GateXYParam_Out.WinSkelPxBlk(1,BlkCt)] = max(SkelPx_Blk_HistCt);
                end
            end
            GateXYParam_Out.WinSkelPxLinIdx_Cell = (WinSkelPxLinIdx_BestMultiWin_Cell{1,SingleWinLengthIdx})'; 
        else
            GateXYParam_Out.WinSkelPxStartIdx = NaN;
            GateXYParam_Out.WinVesselDiameter_px = NaN;
            warning('No Window SkelPxStartIdx found that meets focus requirements.')
        end

    end

end

clearvars SkelPx_Blk_HistCt SkelPx_Blk_HistEdge


%% OPTIONAL: Check that in Single Window Area, Multi Window Case, Windows Do Not Overlap

% % % if (GateXY_ByBlk_Option_In==0) & (GateXY_SingleWinLength_Option_In==1)
% % %     cmap = colormap(jet(numel(GateXYParam_Out.WinSkelPxLinIdx_Cell)+1));
% % %     Img_Temp = zeros(height(ImgStack_In),width(ImgStack_In),'uint8');
% % %     for wCt = 1:numel(GateXYParam_Out.WinSkelPxLinIdx_Cell) % # Windows
% % %         Img_Temp(GateXYParam_Out.WinSkelPxLinIdx_Cell{1,wCt}) = wCt;
% % %     end
% % %     figure;
% % %     imshow(ind2rgb(Img_Temp,cmap),'InitialMagnification',200);
% % % end

clearvars Img_Temp wCt cmap;


%%


clearvars Win_FocusDiameterIdxMean_Max_Temp Win_FocusDiameterIdxMean_MaxIdx_Temp MinIdx_Temp;
clearvars NumSkel_SkelWin_Allwin_Cell;
clearvars SkelWin_Min_tMidPt_Vessel_FlowDiameter_px SkelWin_Skel_bdDistAv_px SkelWin_LinIdx;
clearvars iCt aCt;
