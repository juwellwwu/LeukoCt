clc; close all hidden; fclose('all');
% clearvars;


%% Save list of variables in Workspace after previous module, if does not exist

if ~exist('varList_MOD05') && ~exist('varList_MOD06') % haven't run MOD4 yet
    varList_MOD05 = who;
    varList_MOD05 = varList_MOD05(~ismember(varList_MOD05,varList_All));
    varList_All = vertcat(varList_All,varList_MOD05);
end


%% Extract Reference CBC WBC Concentration (K/uL) for Patient

CellConc_CBC_Tbl = readtable(CBC_FilenameString);

expr_PatientID = 'P\d+';
[StartIdx_PatientID,EndIdx_PatientID] = regexp(SampleIDString,expr_PatientID);

PatientIDString = SampleIDString(StartIdx_PatientID+1:EndIdx_PatientID);
RefConc_KuL = CellConc_CBC_Tbl.RefConc_KuL(CellConc_CBC_Tbl.PatientID==str2num(PatientIDString));

clearvars expr_PatientID StartIdx_PatientID EndIdx_PatientID;


%% Create Output Folder

SaveMODFilePath = strcat(SaveFilePath,'STEP06_WBCConc/');
if exist(SaveMODFilePath,'dir')==7
    rmdir(SaveMODFilePath,'s'); % important
end
mkdir(SaveMODFilePath);


%% Calculate Cell Concentrations (K/uL)

fprintf('Calculating WBC Cell Concentration (K/uL), using the Highest-Z-Score CellCt Window and its Vessel Block) ...\n');

oCellConc_mCell = cell(1, numel(oCellCt.Method));

BrightBlkorWinIdx_mCell = cell(1, numel(oCellCt.Method));
BlkSCS_BrightBlkorWinIdx_mCell = cell(1, numel(oCellCt.Method));

for mCt = 4:4

    if oCellCt.Method(1,mCt)

        % = Specify cell count's window and volume's vessel block
        BrightBlkorWinIdx_mCell{1,mCt} = find(oCellCt.Ct_mCell{1,mCt}.CellCt.BrightRank==1);
        BlkSCS_BrightBlkorWinIdx_mCell{1,mCt} = oCellCt.Ct_mCell{1,mCt}.CellCt.WinSkelPxBlk(BrightBlkorWinIdx_mCell{1,mCt});

        oCellConc = struct;

        % = Incremental cell concentration (K/uL)
        % = i. Volume from Velocity measured using Orient Bands as source of
        % variation (K/uL)
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Cell = cell(1,numel(oCellCt.Method));
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Lo_Cell = cell(1,numel(oCellCt.Method)); 
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Hi_Cell = cell(1,numel(oCellCt.Method));

        % = ii. Volume from Velocity measured using Skeletons as source of
        % variation  (K/uL)
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Cell = cell(1,numel(oCellCt.Method));
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell = cell(1,numel(oCellCt.Method)); 
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell = cell(1,numel(oCellCt.Method)); 

        % = Cumulative cell concentration (K/uL)
        % = i. Volume from Velocity measured using Orient Bands as source of
        % variation (K/uL)
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Cell = cell(1,numel(oCellCt.Method));
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Lo_Cell = cell(1,numel(oCellCt.Method)); 
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Hi_Cell = cell(1,numel(oCellCt.Method)); 

        % = ii. Volume from Velocity measured using Skeletons as source of
        % variation  (K/uL)
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Cell = cell(1,numel(oCellCt.Method));
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell = cell(1,numel(oCellCt.Method));
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell = cell(1,numel(oCellCt.Method));

        % === Incremental Cell Concentration: (# WBC) / (Blood_Vol_per_tMidPt)
        % = i. Volume from Velocity measured using Orient Bands as source of
        % variation (K/uL)
        Vol_AllSkOBwMean = cell2mat(oVolume.V.tMidPt_AllSkOBwMean_Cell);
        Vol_AllSkOBwMean = Vol_AllSkOBwMean(:,:,BlkSCS_BrightBlkorWinIdx_mCell{1,mCt});
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Cell = ... % matrix
            cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}})./Vol_AllSkOBwMean./XY_PxLength.^3.*1E6;
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Cell = ... % convert to cell
            mat2cell(oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Cell,1,repelem(1,NumMidPt),1);

        % Lower concentration limit (larger volume from std error)
        Vol_AllSkOBwMean_Hi = cell2mat(oVolume.V.tMidPt_AllSkOBwMean_Cell) + cell2mat(oVolume.V.tMidPt_AllSkOBwStdEr_Cell); 
        Vol_AllSkOBwMean_Hi = Vol_AllSkOBwMean_Hi(:,:,BlkSCS_BrightBlkorWinIdx_mCell{1,mCt}); 
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Lo_Cell = ...
            cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}})./Vol_AllSkOBwMean_Hi./XY_PxLength.^3.*1E6;
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Lo_Cell = ... % Convert to cell
            mat2cell(oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Lo_Cell,1,repelem(1,NumMidPt),1); 

        % Upper concentration limit (smaller volume from std error)
        Vol_AllSkOBwMean_Lo = cell2mat(oVolume.V.tMidPt_AllSkOBwMean_Cell) - cell2mat(oVolume.V.tMidPt_AllSkOBwStdEr_Cell); 
        Vol_AllSkOBwMean_Lo = Vol_AllSkOBwMean_Lo(:,:,BlkSCS_BrightBlkorWinIdx_mCell{1,mCt}); 
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Hi_Cell = ...
            cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}})./Vol_AllSkOBwMean_Lo./XY_PxLength.^3.*1E6;
        oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Hi_Cell = ...
            mat2cell(oCellConc.ic_KuL.tMidPt_AllSkOBwMean_Hi_Cell,1,repelem(1,NumMidPt),1); 

        % = ii. Volume from Velocity measured using Skeletons as source of
        % variation  (K/uL)
        Vol_1SkOBwMean_SkMean = cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkMean_Cell);
        Vol_1SkOBwMean_SkMean = Vol_1SkOBwMean_SkMean(:,:,BlkSCS_BrightBlkorWinIdx_mCell{1,mCt}); 
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Cell = ...
            cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}})./Vol_1SkOBwMean_SkMean./XY_PxLength.^3.*1E6;
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Cell = ...
            mat2cell(oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Cell,1,repelem(1,NumMidPt),1); 

        % Lower concentration limit (larger volume from std error)
        Vol_1SkOBwMean_SkMean_Hi = cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkMean_Cell) + cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkStdEr_Cell); 
        Vol_1SkOBwMean_SkMean_Hi = Vol_1SkOBwMean_SkMean_Hi(:,:,BlkSCS_BrightBlkorWinIdx_mCell{1,mCt}); 
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell = ...
            cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}})./Vol_1SkOBwMean_SkMean_Hi./XY_PxLength.^3.*1E6;
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell = ...
            mat2cell(oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell,1,repelem(1,NumMidPt),1); 

        % Upper concentration limit (smaller volume from std error)
        Vol_1SkOBwMean_SkMean_Lo = cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkMean_Cell) - cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkStdEr_Cell); % matrix
        Vol_1SkOBwMean_SkMean_Lo = Vol_1SkOBwMean_SkMean_Lo(:,:,BlkSCS_BrightBlkorWinIdx_mCell{1,mCt});  
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell = ...
            cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}})./Vol_1SkOBwMean_SkMean_Lo./XY_PxLength.^3.*1E6;
        oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell = ...
            mat2cell(oCellConc.ic_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell,1,repelem(1,NumMidPt),1); 

        % === Accumulative Cell Concentration: (Accum # WBC) / (Accum_tMidPt_Blood_Vol)
        % = i. Volume from Velocity measured using Orient Bands as source of
        % variation (K/uL)
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Cell = ... % matrix
            cumsum(cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}}))./cumsum(Vol_AllSkOBwMean)./XY_PxLength.^3.*1E6;
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Cell = ... 
            mat2cell(oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Cell,1,repelem(1,NumMidPt),1); 

        % Lower concentration limit (larger volume from std error)
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Lo_Cell = ...
            cumsum(cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}}))./cumsum(Vol_AllSkOBwMean_Hi)./XY_PxLength.^3.*1E6;
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Lo_Cell = ... 
            mat2cell(oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Lo_Cell,1,repelem(1,NumMidPt),1); 

        % Upper concentration limit (smaller volume from std error)
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Hi_Cell = ...
            cumsum(cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}}))./cumsum(Vol_AllSkOBwMean_Lo)./XY_PxLength.^3.*1E6;
        oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Hi_Cell = ...
            mat2cell(oCellConc.cm_KuL.tMidPt_AllSkOBwMean_Hi_Cell,1,repelem(1,NumMidPt),1); 

        % Concentration (K/uL) from Manual CellCt
        if ~isempty(oCellCt.tFr_Manual) % Manual Count info exist
            
            oCellConc.Man_KuL.R_AllSkOBwMean = ...
                oCellCt.tFr_Manual_R_Ct./sum(Vol_AllSkOBwMean,2)./XY_PxLength.^3.*1E6;
            oCellConc.Man_KuL.R_AllSkOBwMean_Lo = ... 
                oCellCt.tFr_Manual_R_Ct./sum(Vol_AllSkOBwMean_Hi,2)./XY_PxLength.^3.*1E6;
            oCellConc.Man_KuL.R_AllSkOBwMean_Hi = ... 
                oCellCt.tFr_Manual_R_Ct./sum(Vol_AllSkOBwMean_Lo,2)./XY_PxLength.^3.*1E6;

            oCellConc.Man_KuL.NR_AllSkOBwMean = ... 
                oCellCt.tFr_Manual_NR_Ct./sum(Vol_AllSkOBwMean,2)./XY_PxLength.^3.*1E6;
            oCellConc.Man_KuL.NR_AllSkOBwMean_Lo = ... 
                oCellCt.tFr_Manual_NR_Ct./sum(Vol_AllSkOBwMean_Hi,2)./XY_PxLength.^3.*1E6;
            oCellConc.Man_KuL.NR_AllSkOBwMean_Hi = ... 
                oCellCt.tFr_Manual_NR_Ct./sum(Vol_AllSkOBwMean_Lo,2)./XY_PxLength.^3.*1E6;

        else 

            oCellConc.Man_KuL.R_AllSkOBwMean = -1;
            oCellConc.Man_KuL.R_AllSkOBwMean_Lo = -1;
            oCellConc.Man_KuL.R_AllSkOBwMean_Hi = -1;

            oCellConc.Man_KuL.NR_AllSkOBwMean = -1;
            oCellConc.Man_KuL.NR_AllSkOBwMean_Lo = -1;
            oCellConc.Man_KuL.NR_AllSkOBwMean_Hi = -1;

        end

        % = ii. Volume from Velocity measured using Skeletons as source of
        % variation  (K/uL)
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Cell = ...
            cumsum(cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}}))./cumsum(Vol_1SkOBwMean_SkMean)./XY_PxLength.^3.*1E6;
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Cell = ...
            mat2cell(oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Cell,1,repelem(1,NumMidPt),1); 

        % Lower concentration limit (larger volume from std error)
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell = ...
            cumsum(cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}}))./cumsum(Vol_1SkOBwMean_SkMean_Hi)./XY_PxLength.^3.*1E6;
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell = ...
            mat2cell(oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Lo_Cell,1,repelem(1,NumMidPt),1); 

        % Upper concentration limit (smaller volume from std error)
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell = ...
            cumsum(cell2mat(oCellCt.Ct_mCell{1,mCt}.CellCt.tMidPt_Ct_Cell{1,BrightBlkorWinIdx_mCell{1,mCt}}))./cumsum(Vol_1SkOBwMean_SkMean_Lo)./XY_PxLength.^3.*1E6;
        oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell = ...
            mat2cell(oCellConc.cm_KuL.tMidPt_1SkOBwMean_SkMean_Hi_Cell,1,repelem(1,NumMidPt),1); 

        % Concentration (K/uL) from Manual Cell Ct 
        if ~isempty(oCellCt.tFr_Manual) % Manual Count info exist
            
            oCellConc.Man_KuL.R_1SkOBwMean_SkMean = ... 
                oCellCt.tFr_Manual_R_Ct./sum(Vol_1SkOBwMean_SkMean,2)./XY_PxLength.^3.*1E6;
            oCellConc.Man_KuL.R_1SkOBwMean_SkMean_Lo = ... 
                oCellCt.tFr_Manual_R_Ct./sum(Vol_1SkOBwMean_SkMean_Hi,2)./XY_PxLength.^3.*1E6;
            oCellConc.Man_KuL.R_1SkOBwMean_SkMean_Hi = ... 
                oCellCt.tFr_Manual_R_Ct./sum(Vol_1SkOBwMean_SkMean_Lo,2)./XY_PxLength.^3.*1E6;

            oCellConc.Man_KuL.NR_1SkOBwMean_SkMean = ... 
                oCellCt.tFr_Manual_NR_Ct./sum(Vol_1SkOBwMean_SkMean,2)./XY_PxLength.^3.*1E6;
            oCellConc.Man_KuL.NR_1SkOBwMean_SkMean_Lo = ... 
                oCellCt.tFr_Manual_NR_Ct./sum(Vol_1SkOBwMean_SkMean_Hi,2)./XY_PxLength.^3.*1E6;
            oCellConc.Man_KuL.NR_1SkOBwMean_SkMean_Hi = ... 
                oCellCt.tFr_Manual_NR_Ct./sum(Vol_1SkOBwMean_SkMean_Lo,2)./XY_PxLength.^3.*1E6;

        else % Fill for reporting

            oCellConc.Man_KuL.R_1SkOBwMean_SkMean = -1;
            oCellConc.Man_KuL.R_1SkOBwMean_SkMean_Lo = -1;
            oCellConc.Man_KuL.R_1SkOBwMean_SkMean_Hi = -1;

            oCellConc.Man_KuL.NR_1SkOBwMean_SkMean = -1;
            oCellConc.Man_KuL.NR_1SkOBwMean_SkMean_Lo = -1;
            oCellConc.Man_KuL.NR_1SkOBwMean_SkMean_Hi = -1;

        end

        % = Prepare output structure 
        oCellConc_mCell{1,mCt} = oCellConc; 
        oCellConc_mCell{1,mCt}.PatientIDStr = PatientIDString;     
        oCellConc_mCell{1,mCt}.Ref_KuL = RefConc_KuL;
        oCellConc_mCell{1,mCt}.BrightBlkorWinIdx = BrightBlkorWinIdx_mCell{1,mCt}; 
        oCellConc_mCell{1,mCt}.BlkSCS_BrightBlkorWinIdx = BlkSCS_BrightBlkorWinIdx_mCell{1,mCt};
        
    end

end

clearvars Vol_* mCt oCellConc;
clearvars PatientIDString RefConc_KuL BrightBlkorWinIdx_mCell BlkSCS_BrightBlkorWinIdx_mCell;


%% Plot Cell Concentrations

for mCt = 4:4

    if oCellCt.Method(1,mCt)

        if mCt==4 
            fprintf('Plotting Cell Concenration (K/uL): METHOD 4) ...\n');
            Method_Str = 'PxWin';
        end

        FigHandle = figure('Position',[20 300 NumMidPt*30+200 375],'visible','on');
        tlo = tiledlayout(2,1,'TileSpacing','Compact','Padding','tight'); 

        cmap = flipud(colormap(abyss(2)));

        % = Subplot 1: incremental cell concentration (#/nL)
        ax = nexttile(tlo);
        hold(ax, 'on');
        grid on;

        x_plot = TimeParam.tSeg_FrMidPt; 
        y_plot = cell2mat(oCellConc_mCell{1,mCt}.ic_KuL.tMidPt_AllSkOBwMean_Cell);
        yl_plot = cell2mat(oCellConc_mCell{1,mCt}.ic_KuL.tMidPt_AllSkOBwMean_Lo_Cell); 
        yh_plot = cell2mat(oCellConc_mCell{1,mCt}.ic_KuL.tMidPt_AllSkOBwMean_Hi_Cell); 

        plot(x_plot,y_plot,'LineWidth',1.5,'Color',cmap(1,:));
        plot(x_plot,yl_plot,'--','LineWidth',1.5,'Color',cmap(1,:));
        plot(x_plot,yh_plot,'--','LineWidth',1.5,'Color',cmap(1,:));

        % Plot Flag Lines
        fdname = fieldnames(tFr_Flag_Plot_mCell{1,mCt});
        FlagID = cellfun(@isempty,strfind(fdname, 'Flag')); 
        fdname(FlagID) = [];
        yf_plot_height = ... 
            linspace(min(yl_plot,[],'all')-0.05*range(y_plot,'all'),min(yl_plot,[],'all')-0.20*range(y_plot,'all'),height(fdname));
        for fgCt = 1:height(fdname) 
            yf_plot = getfield(tFr_Flag_Plot_mCell{1,mCt},fdname{fgCt}).Data;
            cf_plot = getfield(tFr_Flag_Plot_mCell{1,mCt},fdname{fgCt}).Color;
            yf_plot = yf_plot.*yf_plot_height(fgCt);
            yf_plot(abs(yf_plot)<eps) = NaN;
            plot(1:TimeParam.tSeg_FrEnd(end),yf_plot,'LineWidth',2,'Color',cf_plot);
        end

        xlabel('tMidPt_FrCtr (tFr)','Interpreter','none');
        ylabel('Incr WBC Conc (K/uL)','Interpreter','none');
        xlim([0 TimeParam.tSeg_FrEnd(end)]); 
        xticks(0:100:TimeParam.tSeg_FrEnd(end));
        if max(yh_plot,[],'all')>0 
            ylim([0.5*floor(min(yf_plot_height,[],'all')/0.5) 0.5*ceil(max(yh_plot,[],'all')/0.5)]);
        else
            ylim([0.5*floor(min(yf_plot_height,[],'all')/0.5) 0.5]);
        end

        set(gca,'TickDir','out');

        title({strcat('Incremental Cell Conc (K/uL) based on Orient Band Slope Angle Statistics'),...
            strcat('Denominator: oVolume.VR.tMidPt_AllSkOBwMean_Cell +/- oVolume.VR.tMidPt_AllSkOBwStdEr_Cell.'),...
            strcat('CellCt Window @ Vessel Blk:',32,num2str(oCellConc_mCell{1,mCt}.BlkSCS_BrightBlkorWinIdx),':',32,...
            'Flag Lines @ Best CellCt Window & its Vessel Block: Top-to-Bttm: Radius (red), VctyEr.OB (blue), VctyEr.SK (green), FlowStb (purple), RBCBkgd (orange), CellCt (yellow)')},...
            'FontSize',8,'FontWeight','Normal','Interpreter','none');

        % = Subplot 2: cumulative cell concentration (#/nL)
        ax2 = nexttile(tlo);
        hold(ax2, 'on');
        grid on;

        x_plot = TimeParam.tSeg_FrMidPt; 
        y_plot = cell2mat(oCellConc_mCell{1,mCt}.cm_KuL.tMidPt_AllSkOBwMean_Cell);
        yl_plot = cell2mat(oCellConc_mCell{1,mCt}.cm_KuL.tMidPt_AllSkOBwMean_Lo_Cell); 
        yh_plot = cell2mat(oCellConc_mCell{1,mCt}.cm_KuL.tMidPt_AllSkOBwMean_Hi_Cell);

        plot(x_plot,y_plot,'LineWidth',1.5,'Color',cmap(2,:));
        plot(x_plot,yl_plot,'--','LineWidth',1.5,'Color',cmap(2,:));
        plot(x_plot,yh_plot,'--','LineWidth',1.5,'Color',cmap(2,:));

         if ~isempty(oCellCt.tFr_Manual) 
            scatter(TimeParam.tMidPt_FrEnd(end),oCellConc_mCell{1,mCt}.Man_KuL.NR_AllSkOBwMean,10,"red","filled");
            scatter(TimeParam.tMidPt_FrEnd(end),oCellConc_mCell{1,mCt}.Man_KuL.NR_AllSkOBwMean,6,"red","x"');
            scatter(TimeParam.tMidPt_FrEnd(end),oCellConc_mCell{1,mCt}.Man_KuL.NR_AllSkOBwMean,6,"red","x");
        end

        xlabel('tMidPt_FrCtr (tFr)','Interpreter','none');
        ylabel('Cum WBC Conc (K/uL)','Interpreter','none');
        xlim([0 TimeParam.tSeg_FrEnd(end)]); 
        xticks(0:100:TimeParam.tSeg_FrEnd(end));
        if max(yh_plot,[],'all')>0 
            ylim([0.5*floor(min(yl_plot,[],'all')/0.5) 0.5*ceil(max(yh_plot,[],'all')/0.5)]);
        else
            ylim([0.5*floor(min(yl_plot,[],'all')/0.5) 0.5]);
        end
        set(gca,'TickDir','out');

        title({strcat('Cumulative Cell Conc (K/uL) based on Orient Band Slope Angle Statistics. Ref Conc:',32,num2str(oCellConc_mCell{1,mCt}.Ref_KuL,'%0.2f'),'K/uL; Manual Conc (No Repeat):',32,num2str(oCellConc_mCell{1,mCt}.Man_KuL.NR_AllSkOBwMean,'%0.2f'),'K/uL'),...
            strcat('Denominator: (Cumulative) oVolume.VR.tMidPt_AllSkOBwMean_Cell +/- oVolume.VR.tMidPt_AllSkOBwStdEr_Cell.'),...
            strcat('CellCt Window @ Vessel Blk: ',32,num2str(oCellConc_mCell{1,mCt}.BlkSCS_BrightBlkorWinIdx),'; VctyEr Metric (=2*StdErVelocity/MeanVelocity) Mean: ',32,num2str(oFlag_VctyEr.OB.tMidPtwMean_1BlkAllSk_Metric_Cell{1,1,oCellConc_mCell{1,mCt}.BlkSCS_BrightBlkorWinIdx}),',',32,...
            'Median: ',32,num2str(oFlag_VctyEr.OB.tMidPtwMedian_1BlkAllSk_Metric_Cell{1,1,oCellConc_mCell{1,mCt}.BlkSCS_BrightBlkorWinIdx}))},...
            'FontSize',8,'FontWeight','Normal','Interpreter','none');

        % = Save
        exportgraphics(FigHandle,strcat(SaveMODFilePath,Method_Str,'_Plot_tMidPt_WBCConc_K_uL_',timestamp,'.tif'),'Resolution',150);

        close all;
        clearvars FigHandle;

    end

end

clearvars tlo ax* *_plot* cmap fgCt mCt fdname FlagID Method_Str;


%% Clear variables not relevant after this module

clearvars CellConc_CBC_Tbl;


%% Save list of variables in Workspace after previous module, if does not exist

% = OPTION 2: Save only new variables introduced in MOD06
if ~exist('varList_MOD06') && ~exist('varList_MOD07') 
    varList_MOD06 = who;
    varList_MOD06 = varList_MOD06(~ismember(varList_MOD06,varList_All));
    varList_All = vertcat(varList_All,varList_MOD06);
    save(strcat(SaveMODFilePath,'Workspace_MOD06Var_',SampleIDString,'_',timestamp,'.mat'),varList_MOD06{:}); 
end

