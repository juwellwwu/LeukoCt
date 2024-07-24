function [oOrient_Out,oVelocity_Out] = ...
    fCalcOrientVelocity(rp_FFT_Filtered_Clean_Cell_In,oSCS_In,TimeParam_In,Flag_Cell_In)

NumSkel_In = size(rp_FFT_Filtered_Clean_Cell_In,1);
NumSeg_In = size(rp_FFT_Filtered_Clean_Cell_In,2);
NumBlk_In = size(rp_FFT_Filtered_Clean_Cell_In,3);

SkelBlockLength_In = oSCS_In.SkelBlockLength;


%% If "Flag_Cell_In" is empty, create all-zero (no flag) cell

if isempty(Flag_Cell_In)
    Flag_Cell_In = repmat({[0]},NumSkel_In,NumSeg_In,NumBlk_In);
end


%% Calculate Slope Angle and their statistics in Skeleton Coordinate System
% Reference for formula for Weighted Statistics (Case I):
% Bevington, P. R., Data Reduction and Error Analysis for the Physical Sciences, 336 pp., McGraw-Hill, 1969.
% http://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf

fprintf('Calculating Slope Angle and their statistics in Skeleton Coordinate System...\n');

Orient_wMean_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In); 
Orient_wStdev_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In); 
Orient_wStdEr_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In); 
Orient_OBct_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In);

for SegCt = 1:NumSeg_In
    for SkelCt = 1:NumSkel_In
        for BlkCt = 1:NumBlk_In

            if ~(Flag_Cell_In{SkelCt,SegCt,BlkCt})

                rp = rp_FFT_Filtered_Clean_Cell_In{SkelCt,SegCt,BlkCt};
              
                if ~isempty(rp)
                    n_eff = (sum(cat(1,rp.MeanIntensity),1).^2)/sum((cat(1,rp.MeanIntensity)).^2,1); 
                    Orient_wMean_Cell{SkelCt,SegCt,BlkCt} = ...
                        sum(cat(1,rp.Orientation).*cat(1,rp.MeanIntensity),1)./sum(cat(1,rp.MeanIntensity),1);
                    Orient_wStdev_Cell{SkelCt,SegCt,BlkCt} = ...
                        sqrt( (sum(cat(1,rp.MeanIntensity).*(cat(1,rp.Orientation)-Orient_wMean_Cell{SkelCt,SegCt,BlkCt}).^2,1))./(sum(cat(1,rp.MeanIntensity),1)).*(n_eff)./(n_eff-1) ); % 20240201; formula 5
                    Orient_wStdEr_Cell{SkelCt,SegCt,BlkCt} = ...
                        sqrt( (sum(cat(1,rp.MeanIntensity).*(cat(1,rp.Orientation)-Orient_wMean_Cell{SkelCt,SegCt,BlkCt}).^2,1))./(sum(cat(1,rp.MeanIntensity),1)).*(1)./(n_eff-1) ); % 20240201; formula 6
                    Orient_OBct_Cell{SkelCt,SegCt,BlkCt} = height(rp);
                else
                    Orient_wMean_Cell{SkelCt,SegCt,BlkCt} = single(NaN);
                    Orient_wStdev_Cell{SkelCt,SegCt,BlkCt} = single(NaN);
                    Orient_wStdEr_Cell{SkelCt,SegCt,BlkCt} = single(NaN);
                    Orient_OBct_Cell{SkelCt,SegCt,BlkCt} = single(0);
                end

            else

                Orient_wMean_Cell{SkelCt,SegCt,BlkCt} = single(NaN);
                Orient_wStdev_Cell{SkelCt,SegCt,BlkCt} = single(NaN);
                Orient_wStdEr_Cell{SkelCt,SegCt,BlkCt} = single(NaN);
                Orient_OBct_Cell{SkelCt,SegCt,BlkCt} = single(0);

            end

        end
    end
end

clearvars SegCt SkelCt BlkCt rp n_eff;


%% Calculate Flow Velocity (px/fr) and their statistics in Skeleton Coordinate System

fprintf('Calculating Flow Velocity (px/fr) and their statistics in Skeleton Coordinate System...\n');

Velocity_wMean_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In); % in px/fr
Velocity_wStdev_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In); % in px/fr
Velocity_wStdEr_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In); % in px/fr

for SegCt = 1:NumSeg_In
    for SkelCt = 1:NumSkel_In
        for BlkCt = 1:NumBlk_In

            Velocity_wMean_Cell{SkelCt,SegCt,BlkCt} = 1./tand(abs(Orient_wMean_Cell{SkelCt,SegCt,BlkCt}));
            Velocity_wStdev_Cell{SkelCt,SegCt,BlkCt} = 1./tand(abs(Orient_wStdev_Cell{SkelCt,SegCt,BlkCt}));
            Velocity_wStdEr_Cell{SkelCt,SegCt,BlkCt} = 1./tand(abs(Orient_wStdEr_Cell{SkelCt,SegCt,BlkCt}));

        end
    end
end

Velocity_OBct_Cell = Orient_OBct_Cell;

clearvars SegCt SkelCt BlkCt;


%% Calculate tMidPt Values: Assign orientation bands (= slope lines) to appropriate tMidPt

tMidPt_1SkOBCt_Cell = repmat({single(0)},NumSkel_In,TimeParam_In.Num_tMidPt,NumBlk_In); % 1SkOB = Slope lines from 1 skeleton only
tMidPt_1SkOBOrient_Cell = repmat({single([])},NumSkel_In,TimeParam_In.Num_tMidPt,NumBlk_In);
tMidPt_1SkOBWeight_Cell = repmat({single([])},NumSkel_In,TimeParam_In.Num_tMidPt,NumBlk_In);

tMidPt_AllSkOBCt_Cell = repmat({single(0)},1,TimeParam_In.Num_tMidPt,NumBlk_In);  % AllSkOB = Slope lines from all skeletons
tMidPt_AllSkOBOrient_Cell = repmat({single([])},1,TimeParam_In.Num_tMidPt,NumBlk_In);
tMidPt_AllSkOBWeight_Cell = repmat({single([])},1,TimeParam_In.Num_tMidPt,NumBlk_In);

% = Determine max Orientation Band weight
Cell_Temp = reshape(rp_FFT_Filtered_Clean_Cell_In,numel(rp_FFT_Filtered_Clean_Cell_In),1,1);
FlagID = cellfun(@(x) isempty(x),Cell_Temp,'UniformOutput',true);
Cell_Temp(FlagID) = [];
MaxOBMeanIntensity = single(max(cellfun(@(x) max(cat(1,x.MeanIntensity),[],'all'),Cell_Temp,'UniformOutput',true),[],'all'));

clearvars Cell_Temp FlagID;

tMidPt_HistEdge = ... 
    horzcat(TimeParam_In.tMidPt_FrStart,TimeParam_In.tMidPt_FrEnd(end));

for SegCt = 1:NumSeg_In
    for SkelCt = 1:NumSkel_In
        for BlkCt = 1:NumBlk_In

            if ~(Flag_Cell_In{SkelCt,SegCt,BlkCt})

                rp = rp_FFT_Filtered_Clean_Cell_In{SkelCt,SegCt,BlkCt};

                for obCt = 1:height(rp) 
                    if rp(obCt).MeanIntensity/MaxOBMeanIntensity > (1E-6)

                        % = Determine tFr (referenced to full video) and tMidPt occupied by line.
                        [OB_tFr,~] = ind2sub([TimeParam_In.tSeg_FrLength SkelBlockLength_In],rp(obCt).PixelIdxList);
                        OB_tFr = OB_tFr+TimeParam_In.tSeg_FrStart(SegCt)-1; 

                        [~,~,In_tFr_tMidPtBin] = histcounts(OB_tFr,tMidPt_HistEdge);
                        In_tFr_tMidPtBin = unique(In_tFr_tMidPtBin);

                        for binCt = 1:numel(In_tFr_tMidPtBin)

                            % = Collect orient bands from each skeleton separately
                            tMidPt_1SkOBCt_Cell{SkelCt,In_tFr_tMidPtBin(binCt),BlkCt} = ... 
                                tMidPt_1SkOBCt_Cell{SkelCt,In_tFr_tMidPtBin(binCt),BlkCt}+1;
                            tMidPt_1SkOBOrient_Cell{SkelCt,In_tFr_tMidPtBin(binCt),BlkCt} = ...
                                vertcat(tMidPt_1SkOBOrient_Cell{SkelCt,In_tFr_tMidPtBin(binCt),BlkCt},single(rp(obCt).Orientation));
                            tMidPt_1SkOBWeight_Cell{SkelCt,In_tFr_tMidPtBin(binCt),BlkCt} = ...
                                vertcat(tMidPt_1SkOBWeight_Cell{SkelCt,In_tFr_tMidPtBin(binCt),BlkCt},single(rp(obCt).MeanIntensity));

                            % = Collect orient bands from all skeletons
                            tMidPt_AllSkOBCt_Cell{1,In_tFr_tMidPtBin(binCt),BlkCt} = ... 
                                tMidPt_AllSkOBCt_Cell{1,In_tFr_tMidPtBin(binCt),BlkCt}+1;
                            tMidPt_AllSkOBOrient_Cell{1,In_tFr_tMidPtBin(binCt),BlkCt} = ...
                                vertcat(tMidPt_AllSkOBOrient_Cell{1,In_tFr_tMidPtBin(binCt),BlkCt},single(rp(obCt).Orientation));
                            tMidPt_AllSkOBWeight_Cell{1,In_tFr_tMidPtBin(binCt),BlkCt} = ...
                                vertcat(tMidPt_AllSkOBWeight_Cell{1,In_tFr_tMidPtBin(binCt),BlkCt},single(rp(obCt).MeanIntensity));

                        end

                    end
                end

            end

        end
    end
end

% === Methods SK & OB: Calculate flow velocity (px/fr)
tMidPt_1Skel_OBVelocity_Cell = cellfun(@(x) 1./tand(abs(x)), tMidPt_1SkOBOrient_Cell, 'UniformOutput', false); 
tMidPt_AllSkel_OBVelocity_Cell = cellfun(@(x) 1./tand(abs(x)), tMidPt_AllSkOBOrient_Cell, 'UniformOutput', false); 

clearvars rp tMidPt_HistEdge OB_tFr binCt MaxOBMeanIntensity SegCt SkelCt BlkCt obCt;


%% Calculate tMidPt values: Orient Bands from each skeleton treated separately ("1SkOB")

tMidPt_Orient_1SkOBwMean_Cell = ... 
    cellfun(@(wt,or) single((sum(wt.*or,1))/(sum(wt,1))), ...
    tMidPt_1SkOBWeight_Cell,tMidPt_1SkOBOrient_Cell, ...
    'UniformOutput',false);

tMidPt_1SkOBneff_Cell = ... 
    cellfun(@(wt) single((sum(wt,1)).^2./(sum(wt.^2,1))),tMidPt_1SkOBWeight_Cell,'UniformOutput',false);

tMidPt_Orient_1SkOBwStdev_Cell = ... 
    cellfun(@(wt,or,wmor,ct) single(sqrt((sum(wt.*(or-wmor).^2,1))./(sum(wt,1)).*ct./(ct-1))), ...
    tMidPt_1SkOBWeight_Cell, tMidPt_1SkOBOrient_Cell, tMidPt_Orient_1SkOBwMean_Cell, tMidPt_1SkOBneff_Cell, ...
    'UniformOutput',false);

tMidPt_Orient_1SkOBwStdEr_Cell = ... 
    cellfun(@(wt,or,wmor,ct) single(sqrt((sum(wt.*(or-wmor).^2,1))./(sum(wt,1)).*1./(ct-1))), ...
    tMidPt_1SkOBWeight_Cell, tMidPt_1SkOBOrient_Cell, tMidPt_Orient_1SkOBwMean_Cell, tMidPt_1SkOBneff_Cell, ...
    'UniformOutput',false);

tMidPt_Velocity_1SkOBwMean_Cell = ... 
    cellfun(@(wt,vc) single((sum(wt.*vc,1))/(sum(wt,1))), ...
    tMidPt_1SkOBWeight_Cell,tMidPt_1Skel_OBVelocity_Cell, ...
    'UniformOutput',false);

tMidPt_Velocity_1SkOBwStdev_Cell = ... 
    cellfun(@(wt,vc,wmvc,ct) single(sqrt((sum(wt.*(vc-wmvc).^2,1))./(sum(wt,1)).*ct./(ct-1))), ...
    tMidPt_1SkOBWeight_Cell, tMidPt_1Skel_OBVelocity_Cell, tMidPt_Velocity_1SkOBwMean_Cell, tMidPt_1SkOBneff_Cell, ...
    'UniformOutput',false);

tMidPt_Velocity_1SkOBwStdEr_Cell = ... 
    cellfun(@(wt,vc,wmvc,ct) single(sqrt((sum(wt.*(vc-wmvc).^2,1))./(sum(wt,1)).*1./(ct-1))), ...
    tMidPt_1SkOBWeight_Cell, tMidPt_1Skel_OBVelocity_Cell, tMidPt_Velocity_1SkOBwMean_Cell, tMidPt_1SkOBneff_Cell, ...
    'UniformOutput',false);


%% Treatment of cell elements without orient bands: Orient Bands from each skeleton treated separately ("1SkOB")

tMidPt_1Sk_ZeroOBFlag_Mtx = cellfun(@(x) x<2,tMidPt_1SkOBCt_Cell,'UniformOutput',true); 
tMidPt_1Sk_ZeroOBFlag_LinIdx = find(tMidPt_1Sk_ZeroOBFlag_Mtx);

Cell_1Sk_All = {tMidPt_Orient_1SkOBwMean_Cell;...
    tMidPt_Orient_1SkOBwStdev_Cell;...
    tMidPt_Orient_1SkOBwStdEr_Cell;...
    tMidPt_Velocity_1SkOBwMean_Cell;...
    tMidPt_Velocity_1SkOBwStdev_Cell;...
    tMidPt_Velocity_1SkOBwStdEr_Cell};

for cCt = 1:numel(Cell_1Sk_All) 
    
    Mtx_Temp = Cell_1Sk_All{cCt,1}; 
    Mtx_Temp(tMidPt_1Sk_ZeroOBFlag_LinIdx) = {single(-9999)}; 
    Mtx_Temp = cell2mat(Mtx_Temp);
    Mtx_Temp(Mtx_Temp==-9999) = NaN;

    for idxCt = 1:numel(tMidPt_1Sk_ZeroOBFlag_LinIdx)
        [Row,Col,k] = ind2sub(size(Mtx_Temp),tMidPt_1Sk_ZeroOBFlag_LinIdx(idxCt));
        tMidPtAll_Temp = Mtx_Temp(Row,:,k); 
        if (Col>1) && (Col<width(tMidPt_1Sk_ZeroOBFlag_Mtx))
            tMidPtAll_Temp = fillmissing(tMidPtAll_Temp,'linear');
        else
            tMidPtAll_Temp = fillmissing(tMidPtAll_Temp,'nearest');
        end
        Mtx_Temp(Row,:,k) = tMidPtAll_Temp;        
    end

    Cell_1Sk_All{cCt,1} = mat2cell(Mtx_Temp,repelem(1,height(Mtx_Temp)),repelem(1,width(Mtx_Temp)),repelem(1,size(Mtx_Temp,3)));

end

tMidPt_Orient_1SkOBwMean_Cell = Cell_1Sk_All{1,1};
tMidPt_Orient_1SkOBwStdev_Cell = Cell_1Sk_All{2,1};
tMidPt_Orient_1SkOBwStdEr_Cell = Cell_1Sk_All{3,1};
tMidPt_Velocity_1SkOBwMean_Cell = Cell_1Sk_All{4,1};
tMidPt_Velocity_1SkOBwStdev_Cell = Cell_1Sk_All{5,1};
tMidPt_Velocity_1SkOBwStdEr_Cell = Cell_1Sk_All{6,1};

clearvars Cell_1Sk_All tMidPt_1Sk_ZeroOBFlag_Mtx tMidPt_1Sk_ZeroOBFlag_LinIdx Mtx_Temp tMidPtAll_Temp Row Col k cCt idxCt;


%% Calculate tMidPt values: Orient Bands from all skeletons treated together ("AllSkOB")

tMidPt_Orient_AllSkOBwMean_Cell = ... % Weighted average of all slope line orientation in tMidPt
    cellfun(@(wt,or) (sum(wt.*or,1))/(sum(wt,1)), ...
    tMidPt_AllSkOBWeight_Cell,tMidPt_AllSkOBOrient_Cell, ...
    'UniformOutput',false);

tMidPtSkelAll_AllSkOBneff_Cell = ... % = n_eff (formula 4)
    cellfun(@(wt) (sum(wt,1)).^2./(sum(wt.^2,1)),tMidPt_AllSkOBWeight_Cell,'UniformOutput',false);

tMidPt_Orient_AllSkOBwStdev_Cell = ... 
    cellfun(@(wt,or,wmor,ct) sqrt((sum(wt.*(or-wmor).^2,1))./(sum(wt,1)).*ct./(ct-1)), ...
    tMidPt_AllSkOBWeight_Cell, tMidPt_AllSkOBOrient_Cell, tMidPt_Orient_AllSkOBwMean_Cell, tMidPtSkelAll_AllSkOBneff_Cell, ...
    'UniformOutput',false);

tMidPt_Orient_AllSkOBwStdEr_Cell = ... 
    cellfun(@(wt,or,wmor,ct) sqrt((sum(wt.*(or-wmor).^2,1))./(sum(wt,1)).*1./(ct-1)), ...
    tMidPt_AllSkOBWeight_Cell, tMidPt_AllSkOBOrient_Cell, tMidPt_Orient_AllSkOBwMean_Cell, tMidPtSkelAll_AllSkOBneff_Cell, ...
    'UniformOutput',false);

tMidPt_Velocity_AllSkOBwMean_Cell = ... % Weighted average of all slope line velocity (px/fr) in tMidPt
    cellfun(@(wt,vc) (sum(wt.*vc,1))/(sum(wt,1)), ...
    tMidPt_AllSkOBWeight_Cell,tMidPt_AllSkel_OBVelocity_Cell, ...
    'UniformOutput',false);

tMidPt_Velocity_AllSkOBwStdev_Cell = ... % Weighted Standard Deviation of all slope line velocity (px/fr) in tMidPt, always + (formula 5)
    cellfun(@(wt,vc,wmvc,ct) sqrt((sum(wt.*(vc-wmvc).^2,1))./(sum(wt,1)).*ct./(ct-1)), ...
    tMidPt_AllSkOBWeight_Cell, tMidPt_AllSkel_OBVelocity_Cell, tMidPt_Velocity_AllSkOBwMean_Cell, tMidPtSkelAll_AllSkOBneff_Cell, ...
    'UniformOutput',false);

tMidPt_Velocity_AllSkOBwStdEr_Cell = ... % Weighted Standard Error of all slope line velocity (px/fr) in tMidPt, always + (formula 6)
    cellfun(@(wt,vc,wmvc,ct) sqrt((sum(wt.*(vc-wmvc).^2,1))./(sum(wt,1)).*1./(ct-1)), ...
    tMidPt_AllSkOBWeight_Cell, tMidPt_AllSkel_OBVelocity_Cell, tMidPt_Velocity_AllSkOBwMean_Cell, tMidPtSkelAll_AllSkOBneff_Cell, ...
    'UniformOutput',false);


%% Treatment of cell elements without orient bands: Orient Bands from each skeleton treated separately ("AllSkOB")

tMidPt_AllSk_ZeroOBFlag_Mtx = cellfun(@(x) x<2,tMidPt_AllSkOBCt_Cell,'UniformOutput',true); 
tMidPt_AllSk_ZeroOBFlag_LinIdx = find(tMidPt_AllSk_ZeroOBFlag_Mtx);

Cell_AllSk_All = {tMidPt_Orient_AllSkOBwMean_Cell;...
    tMidPt_Orient_AllSkOBwStdev_Cell;...
    tMidPt_Orient_AllSkOBwStdEr_Cell;...
    tMidPt_Velocity_AllSkOBwMean_Cell;...
    tMidPt_Velocity_AllSkOBwStdev_Cell;...
    tMidPt_Velocity_AllSkOBwStdEr_Cell};

for cCt = 1:numel(Cell_AllSk_All) 

    Mtx_Temp = Cell_AllSk_All{cCt,1}; 
    Mtx_Temp(tMidPt_AllSk_ZeroOBFlag_LinIdx) = {single(-9999)}; 
    Mtx_Temp = cell2mat(Mtx_Temp);
    Mtx_Temp(Mtx_Temp==-9999) = NaN;

    for idxCt = 1:numel(tMidPt_AllSk_ZeroOBFlag_LinIdx)        
        [Row,Col,k] = ind2sub(size(Mtx_Temp),tMidPt_AllSk_ZeroOBFlag_LinIdx(idxCt));
        tMidPtAll_Temp = Mtx_Temp(Row,:,k); 
        if (Col>1) && (Col<width(tMidPt_AllSk_ZeroOBFlag_Mtx))
            tMidPtAll_Temp = fillmissing(tMidPtAll_Temp,'linear');
        else
            tMidPtAll_Temp = fillmissing(tMidPtAll_Temp,'nearest');
        end
        Mtx_Temp(Row,:,k) = tMidPtAll_Temp;
    end

     Cell_AllSk_All{cCt,1} = mat2cell(Mtx_Temp,repelem(1,height(Mtx_Temp)),repelem(1,width(Mtx_Temp)),repelem(1,size(Mtx_Temp,3)));

end

tMidPt_Orient_AllSkOBwMean_Cell = Cell_AllSk_All{1,1};
tMidPt_Orient_AllSkOBwStdev_Cell = Cell_AllSk_All{2,1};
tMidPt_Orient_AllSkOBwStdEr_Cell = Cell_AllSk_All{3,1};
tMidPt_Velocity_AllSkOBwMean_Cell = Cell_AllSk_All{4,1};
tMidPt_Velocity_AllSkOBwStdev_Cell = Cell_AllSk_All{5,1};
tMidPt_Velocity_AllSkOBwStdEr_Cell = Cell_AllSk_All{6,1};

clearvars Cell_AllSk_All tMidPt_AllSk_ZeroOBFlag_Mtx tMidPt_AllSk_ZeroOBFlag_LinIdx Mtx_Temp tMidPtAll_Temp Row Col k cCt idxCt;


%% Calculate tMidPt values: Orient Bands from each skeletons treated separately ("1SkOB")

tMidPt_Orient_1SkOBwMean_SkMean_Cell = mat2cell(mean(cell2mat(tMidPt_Orient_1SkOBwMean_Cell),1,'omitmissing'),1,repelem(1,NumSeg_In),repelem(1,NumBlk_In));
tMidPt_Orient_1SkOBwMean_SkStdev_Cell  = mat2cell(std(cell2mat(tMidPt_Orient_1SkOBwMean_Cell),0,1,'omitmissing'),1,repelem(1,NumSeg_In),repelem(1,NumBlk_In));
tMidPt_Orient_1SkOBwMean_SkStdEr_Cell = ... 
    mat2cell(std(cell2mat(tMidPt_Orient_1SkOBwMean_Cell),0,1,'omitmissing')./sqrt(sum(~isnan(cell2mat(tMidPt_Orient_1SkOBwMean_Cell)),1)),1,repelem(1,NumSeg_In),repelem(1,NumBlk_In));

tMidPt_Velocity_1SkOBwMean_SkMean_Cell = mat2cell(mean(cell2mat(tMidPt_Velocity_1SkOBwMean_Cell),1,'omitmissing'),1,repelem(1,NumSeg_In),repelem(1,NumBlk_In));
tMidPt_Velocity_1SkOBwMean_SkStdev_Cell = mat2cell(std(cell2mat(tMidPt_Velocity_1SkOBwMean_Cell),0,1,'omitmissing'),1,repelem(1,NumSeg_In),repelem(1,NumBlk_In));
tMidPt_Velocity_1SkOBwMean_SkStdEr_Cell = ... 
    mat2cell(std(cell2mat(tMidPt_Velocity_1SkOBwMean_Cell),0,1,'omitmissing')./sqrt(sum(~isnan(cell2mat(tMidPt_Velocity_1SkOBwMean_Cell)),1)),1,repelem(1,NumSeg_In),repelem(1,NumBlk_In));

clearvars NumSeg_In NumSkel_In NumBlk_In;


%% Calculate (Weighted) Mean of tMidPt velocity values for reporting 

% CASE 1: Orient Bands from each skeletons treated separately ("1SkOB")
tMidPtwMean_Velocity_1SkOBwMean_Cell = ...
    sum(cell2mat(tMidPt_Velocity_1SkOBwMean_Cell).*repmat(TimeParam_In.tMidPt_FrLength,size(tMidPt_Velocity_1SkOBwMean_Cell,1),1,size(tMidPt_Velocity_1SkOBwMean_Cell,3)),2)./ ...
    sum(repmat(TimeParam_In.tMidPt_FrLength,size(tMidPt_Velocity_1SkOBwMean_Cell,1),1,size(tMidPt_Velocity_1SkOBwMean_Cell,3)),2);

tMidPtwMean_Velocity_1SkOBwStdev_Cell  = ...
    sum(cell2mat(tMidPt_Velocity_1SkOBwStdev_Cell).*repmat(TimeParam_In.tMidPt_FrLength,size(tMidPt_Velocity_1SkOBwStdev_Cell,1),1,size(tMidPt_Velocity_1SkOBwStdev_Cell,3)),2)./ ...
    sum(repmat(TimeParam_In.tMidPt_FrLength,size(tMidPt_Velocity_1SkOBwStdev_Cell,1),1,size(tMidPt_Velocity_1SkOBwStdev_Cell,3)),2);

tMidPtwMean_Velocity_1SkOBwStdEr_Cell  = ...
    sum(cell2mat(tMidPt_Velocity_1SkOBwStdEr_Cell).*repmat(TimeParam_In.tMidPt_FrLength,size(tMidPt_Velocity_1SkOBwStdEr_Cell,1),1,size(tMidPt_Velocity_1SkOBwStdEr_Cell,3)),2)./ ...
    sum(repmat(TimeParam_In.tMidPt_FrLength,size(tMidPt_Velocity_1SkOBwStdEr_Cell,1),1,size(tMidPt_Velocity_1SkOBwStdEr_Cell,3)),2);

tMidPtwMean_Velocity_1SkOBwMean_Cell = mat2cell(tMidPtwMean_Velocity_1SkOBwMean_Cell,repelem(1,size(tMidPt_Velocity_1SkOBwMean_Cell,1)),1,repelem(1,size(tMidPt_Velocity_1SkOBwMean_Cell,3)));
tMidPtwMean_Velocity_1SkOBwStdev_Cell = mat2cell(tMidPtwMean_Velocity_1SkOBwStdev_Cell,repelem(1,size(tMidPt_Velocity_1SkOBwStdev_Cell,1)),1,repelem(1,size(tMidPt_Velocity_1SkOBwStdev_Cell,3)));
tMidPtwMean_Velocity_1SkOBwStdEr_Cell = mat2cell(tMidPtwMean_Velocity_1SkOBwStdEr_Cell,repelem(1,size(tMidPt_Velocity_1SkOBwStdEr_Cell,1)),1,repelem(1,size(tMidPt_Velocity_1SkOBwStdEr_Cell,3)));


% CASE 2: Orient Bands from all skeletons treated together ("AllSkOB")
tMidPtwMean_Velocity_AllSkOBwMean_Cell = ...
        sum(cell2mat(tMidPt_Velocity_AllSkOBwMean_Cell).*repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_AllSkOBwMean_Cell,3)),2)./ ...
        sum(repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_AllSkOBwMean_Cell,3)),2);
tMidPtwMean_Velocity_AllSkOBwStdev_Cell = ...
        sum(cell2mat(tMidPt_Velocity_AllSkOBwStdev_Cell).*repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_AllSkOBwStdev_Cell,3)),2)./ ...
        sum(repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_AllSkOBwStdev_Cell,3)),2);
tMidPtwMean_Velocity_AllSkOBwStdEr_Cell = ...
        sum(cell2mat(tMidPt_Velocity_AllSkOBwStdEr_Cell).*repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_AllSkOBwStdEr_Cell,3)),2)./ ...
        sum(repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_AllSkOBwStdEr_Cell,3)),2);

tMidPtwMean_Velocity_AllSkOBwMean_Cell = mat2cell(tMidPtwMean_Velocity_AllSkOBwMean_Cell,1,1,repelem(1,size(tMidPt_Velocity_AllSkOBwMean_Cell,3)));
tMidPtwMean_Velocity_AllSkOBwStdev_Cell = mat2cell(tMidPtwMean_Velocity_AllSkOBwStdev_Cell,1,1,repelem(1,size(tMidPt_Velocity_AllSkOBwStdev_Cell,3)));
tMidPtwMean_Velocity_AllSkOBwStdEr_Cell = mat2cell(tMidPtwMean_Velocity_AllSkOBwStdEr_Cell,1,1,repelem(1,size(tMidPt_Velocity_AllSkOBwStdEr_Cell,3)));


% CASE 3: Orient Bands from each skeletons treated separately ("1SkOB")
tMidPtwMean_Velocity_1SkOBwMean_SkMean_Cell = ...
        sum(cell2mat(tMidPt_Velocity_1SkOBwMean_SkMean_Cell).*repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_1SkOBwMean_SkMean_Cell,3)),2)./ ...
        sum(repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_1SkOBwMean_SkMean_Cell,3)),2);
tMidPtwMean_Velocity_1SkOBwMean_SkStdev_Cell = ...
        sum(cell2mat(tMidPt_Velocity_1SkOBwMean_SkStdev_Cell).*repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_1SkOBwMean_SkStdev_Cell,3)),2)./ ...
        sum(repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_1SkOBwMean_SkStdev_Cell,3)),2);
tMidPtwMean_Velocity_1SkOBwMean_SkStdEr_Cell = ...
        sum(cell2mat(tMidPt_Velocity_1SkOBwMean_SkStdEr_Cell).*repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_1SkOBwMean_SkStdEr_Cell,3)),2)./ ...
        sum(repmat(TimeParam_In.tMidPt_FrLength,1,1,size(tMidPt_Velocity_1SkOBwMean_SkStdEr_Cell,3)),2);

tMidPtwMean_Velocity_1SkOBwMean_SkMean_Cell = mat2cell(tMidPtwMean_Velocity_1SkOBwMean_SkMean_Cell,1,1,repelem(1,size(tMidPt_Velocity_1SkOBwMean_SkMean_Cell,3)));
tMidPtwMean_Velocity_1SkOBwMean_SkStdev_Cell = mat2cell(tMidPtwMean_Velocity_1SkOBwMean_SkStdev_Cell,1,1,repelem(1,size(tMidPt_Velocity_1SkOBwMean_SkStdev_Cell,3)));
tMidPtwMean_Velocity_1SkOBwMean_SkStdEr_Cell = mat2cell(tMidPtwMean_Velocity_1SkOBwMean_SkStdEr_Cell,1,1,repelem(1,size(tMidPt_Velocity_1SkOBwMean_SkStdEr_Cell,3)));


%% Load Output Structure

oOrient_Out = struct;
oVelocity_Out = struct;

% === Orient_Param
% i. SCS: Orient bands from each {SkelCt,SegCt,BlkCt} treated separately 
oOrient_Out.wMean_Cell = Orient_wMean_Cell;
oOrient_Out.wStdev_Cell = Orient_wStdev_Cell;
oOrient_Out.wStdEr_Cell = Orient_wStdEr_Cell;
oOrient_Out.OBct_Cell = Orient_OBct_Cell;

% ii. tMidPt: Orient bands from each skeleton treated separately
oOrient_Out.tMidPt_1SkOBwMean_Cell = tMidPt_Orient_1SkOBwMean_Cell; % size = NumSkel x NumMidPt x NumBlk
oOrient_Out.tMidPt_1SkOBwStdev_Cell = tMidPt_Orient_1SkOBwStdev_Cell ; % size = NumSkel x NumMidPt x NumBlk
oOrient_Out.tMidPt_1SkOBwStdEr_Cell = tMidPt_Orient_1SkOBwStdEr_Cell ; % size = NumSkel x NumMidPt x NumBlk
oOrient_Out.tMidPt_1SkOBct_Cell = tMidPt_1SkOBCt_Cell; % size = NumSkel x NumMidPt x NumBlk

% iii. tMidPt: Orient bands from all skeleton treated together
% For Flow Velocity Measurement: statistics from Orient Band slope variation
oOrient_Out.tMidPt_AllSkOBwMean_Cell = tMidPt_Orient_AllSkOBwMean_Cell; % size = 1 x NumMidPt x NumBlk
oOrient_Out.tMidPt_AllSkOBwStdev_Cell = tMidPt_Orient_AllSkOBwStdev_Cell ; % size = 1 x NumMidPt x NumBlk
oOrient_Out.tMidPt_AllSkOBwStdEr_Cell = tMidPt_Orient_AllSkOBwStdEr_Cell ; % size = 1 x NumMidPt x NumBlk
oOrient_Out.tMidPt_AllSkOBct_Cell = tMidPt_AllSkOBCt_Cell; % size = 1 x NumMidPt x NumBlk

% iv. tMidPt: Orient bands from each skeleton treated separately
oOrient_Out.tMidPt_1SkOBwMean_SkMean_Cell = tMidPt_Orient_1SkOBwMean_SkMean_Cell; % size = 1 x NumMidPt x NumBlk
oOrient_Out.tMidPt_1SkOBwMean_SkStdev_Cell = tMidPt_Orient_1SkOBwMean_SkStdev_Cell; % size = 1 x NumMidPt x NumBlk
oOrient_Out.tMidPt_1SkOBwMean_SkStdEr_Cell = tMidPt_Orient_1SkOBwMean_SkStdEr_Cell; % size = 1 x NumMidPt x NumBlk

% === Velocity_Param
% i. SCS: Orient bands from each {SkelCt,SegCt,BlkCt} treated separately 
oVelocity_Out.wMean_Cell = Velocity_wMean_Cell;
oVelocity_Out.wStdev_Cell = Velocity_wStdev_Cell;
oVelocity_Out.wStdEr_Cell = Velocity_wStdEr_Cell;
oVelocity_Out.OBct_Cell = Velocity_OBct_Cell;

% ii. tMidPt: Orient bands from each skeleton treated separately
oVelocity_Out.tMidPt_1SkOBwMean_Cell = tMidPt_Velocity_1SkOBwMean_Cell; % size = NumSkel x NumMidPt x NumBlk
oVelocity_Out.tMidPt_1SkOBwStdev_Cell = tMidPt_Velocity_1SkOBwStdev_Cell ; % size = NumSkel x NumMidPt x NumBlk
oVelocity_Out.tMidPt_1SkOBwStdEr_Cell = tMidPt_Velocity_1SkOBwStdEr_Cell ; % size = NumSkel x NumMidPt x NumBlk
oVelocity_Out.tMidPt_1SkOBct_Cell = tMidPt_1SkOBCt_Cell; % size = NumSkel x NumMidPt x NumBlk

% (Weighted) mean of tMidPt values for reporting
oVelocity_Out.tMidPtwMean_1SkOBwMean_Cell = tMidPtwMean_Velocity_1SkOBwMean_Cell; % size = NumSkel x 1 x NumBlk
oVelocity_Out.tMidPtwMean_1SkOBwStdev_Cell = tMidPtwMean_Velocity_1SkOBwStdev_Cell; % size = NumSkel x 1 x NumBlk
oVelocity_Out.tMidPtwMean_1SkOBwStdEr_Cell = tMidPtwMean_Velocity_1SkOBwStdEr_Cell; % size = NumSkel x 1 x NumBlk

% iii. tMidPt: Orient bands from all skeleton treated together
% For Flow Velocity Measurement: statistics from Orient Band slope variation
oVelocity_Out.tMidPt_AllSkOBwMean_Cell = tMidPt_Velocity_AllSkOBwMean_Cell; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPt_AllSkOBwStdev_Cell = tMidPt_Velocity_AllSkOBwStdev_Cell ; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPt_AllSkOBwStdEr_Cell = tMidPt_Velocity_AllSkOBwStdEr_Cell ; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPt_AllSkOBct_Cell = tMidPt_AllSkOBCt_Cell; % size = 1 x NumMidPt x NumBlk

% (Weighted) mean of tMidPt values, for reporting
oVelocity_Out.tMidPtwMean_AllSkOBwMean_Cell = tMidPtwMean_Velocity_AllSkOBwMean_Cell; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPtwMean_AllSkOBwStdev_Cell = tMidPtwMean_Velocity_AllSkOBwStdev_Cell ; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPtwMean_AllSkOBwStdEr_Cell = tMidPtwMean_Velocity_AllSkOBwStdEr_Cell ; % size = 1 x NumMidPt x NumBlk

% iv. tMidPt: Orient bands from each skeleton treated separately
oVelocity_Out.tMidPt_1SkOBwMean_SkMean_Cell = tMidPt_Velocity_1SkOBwMean_SkMean_Cell; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPt_1SkOBwMean_SkStdev_Cell = tMidPt_Velocity_1SkOBwMean_SkStdev_Cell; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPt_1SkOBwMean_SkStdEr_Cell = tMidPt_Velocity_1SkOBwMean_SkStdEr_Cell; % size = 1 x NumMidPt x NumBlk

% (Weighted) mean of tMidPt values, for reporting
oVelocity_Out.tMidPtwMean_1SkOBwMean_SkMean_Cell = tMidPtwMean_Velocity_1SkOBwMean_SkMean_Cell; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPtwMean_1SkOBwMean_SkStdev_Cell = tMidPtwMean_Velocity_1SkOBwMean_SkStdev_Cell; % size = 1 x NumMidPt x NumBlk
oVelocity_Out.tMidPtwMean_1SkOBwMean_SkStdEr_Cell = tMidPtwMean_Velocity_1SkOBwMean_SkStdEr_Cell; % size = 1 x NumMidPt x NumBlk



