clc; close all hidden; fclose('all');
% clearvars;


%% Save list of variables in Workspace after previous module, if does not exist

if ~exist('varList_MOD02') && ~exist('varList_MOD04') 
    varList_MOD02 = who;
    varList_MOD02 = varList_MOD02(~ismember(varList_MOD02,varList_All));
     varList_All = vertcat(varList_All,varList_MOD02);
end

if ~exist('DEMOOption')
    DEMOOption = 'n';
end


%% Create Output Folder

SaveMODFilePath = strcat(SaveFilePath,'STEP04_Velocity_Volume/');
if exist(SaveMODFilePath,'dir')==7
    rmdir(SaveMODFilePath,'s'); % important
end
mkdir(SaveMODFilePath);


%% Calculate tMidPt VolRate (px^3/fr) and Volume (px^3)

TimeParam.tMidPt_Cell = ... 
    repmat(mat2cell(single(TimeParam.tMidPt_FrLength),1,repelem(1,NumMidPt)),NumSkel_VsWidth,1,NumBlk);
TimeParam.tMidPt_AllSk_Cell = ... 
    repmat(mat2cell(single(TimeParam.tMidPt_FrLength),1,repelem(1,NumMidPt)),1,1,NumBlk);

% = Calculate tMidPt Flow Volume Rate (px^3/fr)

oVolume = struct;

% i. Statistics with orient band slopes as sample size and source of variation.
oVolume.VR.tMidPt_AllSkOBwMean_Cell = ... 
    cellfun(@times,oVelocity.tMidPt_AllSkOBwMean_Cell,oRadius.L.tMidPt_AllSk_Area_Cell, 'UniformOutput',false);
oVolume.VR.tMidPt_AllSkOBwStdEr_Cell = ... 
    cellfun(@times,oVelocity.tMidPt_AllSkOBwStdEr_Cell,oRadius.L.tMidPt_AllSk_Area_Cell, 'UniformOutput',false);
oVolume.VR.tMidPt_AllSkOBct_Cell = oVelocity.tMidPt_AllSkOBct_Cell; 

% ii. Statistics with skeletons as sample size and source of variation
oVolume.VR.tMidPt_1SkOBwMean_SkMean_Cell = ... 
    cellfun(@times,oVelocity.tMidPt_1SkOBwMean_SkMean_Cell,oRadius.L.tMidPt_AllSk_Area_Cell, 'UniformOutput',false);
oVolume.VR.tMidPt_1SkOBwMean_SkStdEr_Cell = ... 
    cellfun(@times,oVelocity.tMidPt_1SkOBwMean_SkStdEr_Cell,oRadius.L.tMidPt_AllSk_Area_Cell, 'UniformOutput',false);

% = Calculate tMidPt Flow Volume (px^3)
% i. Statistics with orient band slopes as sample size and source of variation.
oVolume.V.tMidPt_AllSkOBwMean_Cell = ... 
    cellfun(@times,oVolume.VR.tMidPt_AllSkOBwMean_Cell,TimeParam.tMidPt_AllSk_Cell, 'UniformOutput',false);
oVolume.V.tMidPt_AllSkOBwStdEr_Cell = ...
    cellfun(@times,oVolume.VR.tMidPt_AllSkOBwStdEr_Cell,TimeParam.tMidPt_AllSk_Cell, 'UniformOutput',false);

% ii. Statistics with skeletons as sample size and source of variation
oVolume.V.tMidPt_1SkOBwMean_SkMean_Cell = ...
    cellfun(@times,oVolume.VR.tMidPt_1SkOBwMean_SkMean_Cell,TimeParam.tMidPt_AllSk_Cell, 'UniformOutput',false);
oVolume.V.tMidPt_1SkOBwMean_SkStdEr_Cell = ... 
    cellfun(@times,oVolume.VR.tMidPt_1SkOBwMean_SkStdEr_Cell,TimeParam.tMidPt_AllSk_Cell, 'UniformOutput',false);

% = Total flow volume for each vessel blk

% i. Statistics with orient band slopes as sample size and source of variation. 
oVolume.V.tMidPtSum_AllSkOBwMean = sum(cell2mat(oVolume.V.tMidPt_AllSkOBwMean_Cell),2); 
oVolume.V.tMidPtSum_AllSkOBwStdEr = sum(cell2mat(oVolume.V.tMidPt_AllSkOBwStdEr_Cell),2);

oVolume.V.tMidPtSum_AllSkOBwMean_nL = oVolume.V.tMidPtSum_AllSkOBwMean.*XY_PxLength.^3./1E6; 
oVolume.V.tMidPtSum_AllSkOBwStdEr_nL = oVolume.V.tMidPtSum_AllSkOBwStdEr.*XY_PxLength.^3./1E6;

% ii. Statistics with skeletons as sample size and source of variation
oVolume.V.tMidPtSum_1SkOBwMean_SkMean = sum(cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkMean_Cell),2);  
oVolume.V.tMidPtSum_1SkOBwMean_SkStdEr = sum(cell2mat(oVolume.V.tMidPt_1SkOBwMean_SkStdEr_Cell),2);

oVolume.V.tMidPtSum_1SkOBwMean_SkMean_nL = oVolume.V.tMidPtSum_1SkOBwMean_SkMean.*XY_PxLength.^3./1E6; 
oVolume.V.tMidPtSum_1SkOBwMean_SkStdEr_nL = oVolume.V.tMidPtSum_1SkOBwMean_SkStdEr.*XY_PxLength.^3./1E6;


%% Locate & Flag time frame stretches with large velocity standard errors

[oFlag_VctyEr] = fCreateVelocityStdErFlagCell(oVelocity,TimeParam,oSCS,0.06,0.10,find(oFocusQ.BlkRank==1));


%% Locate & flag time frame stretches with unstable flow [Removed]

[oFlag_FlowStb] = fCreateDUMMYFlowStabilityFlagCell(oOrient,TimeParam,-9999,oSCS,find(oFocusQ.BlkRank==1)); % Dummy


%% Prepare Flag Struct for Plotting

tFr_Flag_Prelim_Plot = struct;

tFr_Flag_Prelim_Plot.cmap = colormap(slanCM('Pastel1',9)); 
close all; 

tFr_Flag_Prelim_Plot.Flag01.Data = oFlag_Radius.L.tFr_1BlkAllSk_Flag_Cell{1,1,find(oFocusQ.BlkRank==1)}; 
tFr_Flag_Prelim_Plot.Flag01.Color = tFr_Flag_Prelim_Plot.cmap(1,:);

tFr_Flag_Prelim_Plot.Flag02.Data = oFlag_VctyEr.OB.tFr_1BlkAllSk_Flag_Cell{1,1,find(oFocusQ.BlkRank==1)}; 
tFr_Flag_Prelim_Plot.Flag02.Color = tFr_Flag_Prelim_Plot.cmap(2,:);

tFr_Flag_Prelim_Plot.Flag03.Data = oFlag_VctyEr.SK.tFr_1BlkAllSk_Flag_Cell{1,1,find(oFocusQ.BlkRank==1)}; 
tFr_Flag_Prelim_Plot.Flag03.Color = tFr_Flag_Prelim_Plot.cmap(3,:);

tFr_Flag_Prelim_Plot.Flag04.Data = oFlag_FlowStb.tFr_1BlkAllSk_Flag_Cell{1,1,find(oFocusQ.BlkRank==1)}; 
tFr_Flag_Prelim_Plot.Flag04.Color = tFr_Flag_Prelim_Plot.cmap(4,:);


%% Plot tMidPt Flow Velocity (um/s) 

FigHandle1 = figure('Position',[20 300 NumMidPt*30+200 600],'visible','on');
tlo = tiledlayout(3,1,'TileSpacing','Compact','Padding','tight'); 

% = Subplot 1
ax = nexttile(tlo); 
hold(ax, 'on');
grid on;

cmap = colormap(summer(round(NumBlk*1.5))); % Color vessel blocks from dark green to grass green

for BlkCt = 1:NumBlk 

      errorbar(TimeParam.tSeg_FrMidPt,...
        cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell(:,:,BlkCt)).*XY_PxLength./Z_PxLength,...
        cell2mat(oVelocity.tMidPt_AllSkOBwStdEr_Cell(:,:,BlkCt)).*XY_PxLength./Z_PxLength,... 
        cell2mat(oVelocity.tMidPt_AllSkOBwStdEr_Cell(:,:,BlkCt)).*XY_PxLength./Z_PxLength,... 
        'LineWidth',1,'Color',cmap(BlkCt,:));

end

xlabel('Time Frame','Interpreter','none');
ylabel('Velocity (um/s)');

xlim([0 TimeParam.tSeg_FrEnd(end)]); 
xticks(0:100:TimeParam.tSeg_FrEnd(end));
ymin = min((cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell)-cell2mat(oVelocity.tMidPt_AllSkOBwStdEr_Cell)).*XY_PxLength./Z_PxLength,[],'all');
ymax = max((cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell)+cell2mat(oVelocity.tMidPt_AllSkOBwStdEr_Cell)).*XY_PxLength./Z_PxLength,[],'all');
ylim([100*floor(ymin/100) 100*ceil(ymax/100)]);

title({strcat('Blood Flow Velocity (um/s) based on Orient Band Slope Angle Statistics'),...
    strcat('oVelocity.tMidPt_AllSkOBwMean_Cell +/- oVelocity.tMidPt_AllSkOBwStdEr_Cell'),...
    strcat('Best FocusQ Blk: ',32,num2str(find(oFocusQ.BlkRank==1)),'; VctyEr Metric (=2*StdErVelocity/MeanVelocity) Mean: ',32,num2str(oFlag_VctyEr.OB.tMidPtwMean_1BlkAllSk_Metric_Cell{1,1,find(oFocusQ.BlkRank==1)}),',',32,...
    'Median: ',32,num2str(oFlag_VctyEr.OB.tMidPtwMedian_1BlkAllSk_Metric_Cell{1,1,find(oFocusQ.BlkRank==1)})),...
    strcat('Mean Per-Blk Velocity (um/s):',32,char(strjoin(string(mean(cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell).*XY_PxLength./Z_PxLength,2)))),';',32,...
    'Mean All-Blk Velocity:',32,num2str(mean(cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell).*XY_PxLength./Z_PxLength,'all')))},...
    'FontSize',8,'FontWeight','Normal','Interpreter','none');

hold off;

clearvars cmap ymin ymax BlkCt;

% = Subplot 2
ax2 = nexttile(tlo); 
hold(ax2, 'on');
grid on;

cmap = colormap(summer(round(NumBlk*1.5))); 

plot(TimeParam.tSeg_FrMidPt,cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell(:,:,find(oFocusQ.BlkRank==1))).*XY_PxLength./Z_PxLength,...
    ':','LineWidth',1,'Color',cmap(find(oFocusQ.BlkRank==1),:));

fdname = fieldnames(tFr_Flag_Prelim_Plot);
FlagID = cellfun(@isempty,strfind(fdname, 'Flag')); 
fdname(FlagID) = [];
yh = ... 
    linspace(min(cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell(:,:,find(oFocusQ.BlkRank==1))).*XY_PxLength./Z_PxLength,[],'all')-50,...
    min(cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell(:,:,find(oFocusQ.BlkRank==1))).*XY_PxLength./Z_PxLength,[],'all')-200,...
    height(fdname)); 
for fgCt = 1:height(fdname) % # Flags
    y = getfield(tFr_Flag_Prelim_Plot,fdname{fgCt}).Data;
    c = getfield(tFr_Flag_Prelim_Plot,fdname{fgCt}).Color;
    y = y.*yh(fgCt);
    y(y<eps) = NaN;
    plot(1:TimeParam.tSeg_FrEnd(end),y,'LineWidth',2,'Color',c);
end

xlabel('Time Frame','Interpreter','none');
ylabel('Velocity (um/s)');

xlim([0 TimeParam.tSeg_FrEnd(end)]); 
xticks(0:100:TimeParam.tSeg_FrEnd(end));
ymin = min(yh,[],"all");
ymax = max(cell2mat(oVelocity.tMidPt_AllSkOBwMean_Cell(:,:,find(oFocusQ.BlkRank==1))).*XY_PxLength./Z_PxLength,[],'all');
ylim([100*floor(ymin/100) 100*ceil(ymax/100)]);

title({strcat('Flag Lines based on Orient Band Slope Angle Statistics: Best FocusQ Blk ONLY'),...
    strcat('Dotted Line: oVelocity.tMidPt_AllSkOBwMean_Cell'),...
    strcat('Solid Flag Lines: Top-to-Bttm: Radius (red), VctyEr.OB (blue), VctyEr.SK (green), FlowStb (purple)')},...
    'FontSize',8,'FontWeight','Normal','Interpreter','none');

hold off;

clearvars cmap fdname FlagID y y_AllBlk yh yCt c c_AllBlk fgCt fBlkCt ymin ymax BlkCt;

% = Subplot 3
ax3 = nexttile(tlo); 
hold(ax3, 'on');
grid on;

cmap = colormap(summer(round(NumBlk*1.5))); 

for BlkCt = 1:NumBlk 
    scatter(TimeParam.tSeg_FrMidPt,cell2mat(oVelocity.tMidPt_AllSkOBct_Cell(:,:,BlkCt)),30,cmap(BlkCt,:));
end

xlabel('Time Frame','Interpreter','none');
ylabel('OB Ct');

xlim([0 TimeParam.tSeg_FrEnd(end)]); 
xticks(0:100:TimeParam.tSeg_FrEnd(end));
ymin = min(cell2mat(oVelocity.tMidPt_AllSkOBct_Cell),[],'all');
ymax = max(cell2mat(oVelocity.tMidPt_AllSkOBct_Cell),[],'all');
ylim([20*floor(ymin/20) 20*ceil(ymax/20)]);

title({strcat('Orient Band Count contributing to Blood Flow Velocity (um/s) based on Orient Band Slope Angle Statistics')},...
    'FontSize',8,'FontWeight','Normal','Interpreter','none');

hold off;

clearvars cmap ymin ymax BlkCt;

% = Save Plot
exportgraphics(FigHandle1,strcat(SaveMODFilePath,'Plot_tMidPt_Velocity_um_s_',timestamp,'.tif'),'Resolution',150);

close all;
clearvars ax* tlo FigHandle1;


%% Plot tMidPt Flow Volume Rate (nL/min) and Flow Volume (nL) 

FigHandle2 = figure('Position',[70 250 NumMidPt*30+200 400],'visible','on');
tlo = tiledlayout(2,1,'TileSpacing','Compact','Padding','tight'); % tile layout; replace subplot()

% = Subplot 1
ax = nexttile(tlo); 
hold(ax, 'on');
grid on;

cmap = colormap(summer(round(NumBlk*1.5))); 

for BlkCt = 1:NumBlk 

      errorbar(TimeParam.tSeg_FrMidPt,...
        cell2mat(oVolume.VR.tMidPt_AllSkOBwMean_Cell(:,:,BlkCt)).*XY_PxLength.^3./Z_PxLength./1E6.*60,...
        cell2mat(oVolume.VR.tMidPt_AllSkOBwStdEr_Cell(:,:,BlkCt)).*XY_PxLength.^3./Z_PxLength./1E6.*60,... 
        cell2mat(oVolume.VR.tMidPt_AllSkOBwStdEr_Cell(:,:,BlkCt)).*XY_PxLength.^3./Z_PxLength./1E6.*60,... 
        'LineWidth',1,'Color',cmap(BlkCt,:));

end

xlabel('Time Frame','Interpreter','none');
ylabel('Flow Rate (nL/min)');

xlim([0 TimeParam.tSeg_FrEnd(end)]); 
xticks(0:100:TimeParam.tSeg_FrEnd(end));
ymin = min((cell2mat(oVolume.VR.tMidPt_AllSkOBwMean_Cell)-cell2mat(oVolume.VR.tMidPt_AllSkOBwStdEr_Cell)).*XY_PxLength.^3./Z_PxLength./1E6.*60,[],'all');
ymax = max((cell2mat(oVolume.VR.tMidPt_AllSkOBwMean_Cell)+cell2mat(oVolume.VR.tMidPt_AllSkOBwStdEr_Cell)).*XY_PxLength.^3./Z_PxLength./1E6.*60,[],'all');
ylim([5*floor(ymin/5) 5*ceil(ymax/5)]);

title({strcat('Blood Flow Vol Rate (nL/min) based on Orient Band Slope Angle Statistics'),...
    strcat('oVolume.VR.tMidPt_AllSkOBwMean_Cell +/- oVolume.VR.tMidPt_AllSkOBwStdEr_Cell; Best FocusQ Blk = ',32,num2str( find(oFocusQ.BlkRank==1))),...
    strcat('Mean Per-Blk Vol Rate (nL/min):',32,char(strjoin(string(mean(cell2mat(oVolume.VR.tMidPt_AllSkOBwMean_Cell).*XY_PxLength.^3./Z_PxLength./1E6.*60,2)))),';',32,...
    'Mean All-Blk Vol Rate:',32,num2str(mean(cell2mat(oVolume.VR.tMidPt_AllSkOBwMean_Cell).*XY_PxLength.^3./Z_PxLength./1E6.*60,'all')))},...
    'FontSize',8,'FontWeight','Normal','Interpreter','none');

hold off;

clearvars cmap ymin ymax BlkCt;

% = Subplot 2
ax2 = nexttile(tlo); 
hold(ax2, 'on');
grid on;

cmap = colormap(summer(round(NumBlk*1.5)));

for BlkCt = 1:NumBlk 

      errorbar(TimeParam.tSeg_FrMidPt,...
        cell2mat(oVolume.V.tMidPt_AllSkOBwMean_Cell(:,:,BlkCt)).*XY_PxLength.^3./1E6,...
        cell2mat(oVolume.V.tMidPt_AllSkOBwStdEr_Cell(:,:,BlkCt)).*XY_PxLength.^3./1E6,... 
        cell2mat(oVolume.V.tMidPt_AllSkOBwStdEr_Cell(:,:,BlkCt)).*XY_PxLength.^3./1E6,... 
        'LineWidth',1,'Color',cmap(BlkCt,:));

end

xlabel('Time Frame','Interpreter','none');
ylabel('Flow Vol (nL)');

xlim([0 TimeParam.tSeg_FrEnd(end)]); 
xticks(0:100:TimeParam.tSeg_FrEnd(end));
ymin = min((cell2mat(oVolume.V.tMidPt_AllSkOBwMean_Cell)-cell2mat(oVolume.V.tMidPt_AllSkOBwStdEr_Cell)).*XY_PxLength.^3./1E6,[],'all');
ymax = max((cell2mat(oVolume.V.tMidPt_AllSkOBwMean_Cell)+cell2mat(oVolume.V.tMidPt_AllSkOBwStdEr_Cell)).*XY_PxLength.^3./1E6,[],'all');
ylim([0.1*floor(ymin/0.1) 0.1*ceil(ymax/0.1)]);

title({strcat('Blood Flow Volume (nL) based on Orient Band Slope Angle Statistics'),...
    strcat('oVolume.V.tMidPt_AllSkOBwMean_Cell +/- oVolume.V.tMidPt_AllSkOBwStdEr_Cell.'),...
    strcat('Per-Blk Total Vol (nL):',32,char(strjoin(string(oVolume.V.tMidPtSum_AllSkOBwMean_nL))),',',32,...
    'StdEr:',32,char(strjoin(string(oVolume.V.tMidPtSum_AllSkOBwStdEr_nL))),';',32,...
    'Mean Per-Blk Total Vol:',32,num2str(mean(oVolume.V.tMidPtSum_AllSkOBwMean_nL,'all')),',',32,...
    'Mean StdEr:',32,num2str(mean(oVolume.V.tMidPtSum_AllSkOBwStdEr_nL,'all')))},...
    'FontSize',8,'FontWeight','Normal','Interpreter','none');

hold off;

clearvars cmap ymin ymax BlkCt;

% = Save Plot
exportgraphics(FigHandle2,strcat(SaveMODFilePath,'Plot_tMidPt_Volume_nL_min_',timestamp,'.tif'),'Resolution',150);

close all;
clearvars ax* tlo FigHandle2;


%% Clear variables likely no longer relevant after this module

clearvars ImgDR_FFT_Filtered_Clean*;
clearvars rp_FFT_Filtered_Clean*;
clearvars LamFlow*;
clearvars FrMidPtSkelPxAv_Vessel_Area_Cell;
clearvars tFr_RegQuality_Flag tFr_RegQualityMetric *RegQuality*;


%% Save list of variables in Workspace after previous module, if does not exist

% = OPTION 2: Save only new variables introduced
if ~exist('varList_MOD04') && ~exist('varList_MOD05') 
    varList_MOD04 = who;
    varList_MOD04 = varList_MOD04(~ismember(varList_MOD04,varList_All));
    varList_All = vertcat(varList_All,varList_MOD04);
    save(strcat(SaveMODFilePath,'Workspace_STEP04Var_',SampleIDString,'_',timestamp,'.mat'),varList_MOD04{:},...
        'TimeParam','oOrient','oVelocity'); % updated variables
end


%% Save script in directory

ScriptName=mfilename;
PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveMODFilePath);
publish(strcat(ScriptName,'.m'),PublishOptions);

