function [Fr_ManualCellCt_Out] = fLoadManualCellCtFr(Fr_ManualCellCt_FoldernameString_In,SampleIDString_In,ZCropParam_In)

%%  Load Manual Cell Ct Frame information, if info exists in directory

fprintf('Loading Manual Cell Ct Frame information, if present ...\n');

Fr_ManualCellCt_FilenameString = strcat(Fr_ManualCellCt_FoldernameString_In,SampleIDString_In,'_ManualCellCt');

% = Read frames of cell appearance for manual cell count
if isfile(Fr_ManualCellCt_FilenameString)    
    fileID_Fr_ManualCellCt = fopen(Fr_ManualCellCt_FilenameString,'r');

    formatSpec_Fr_ManualCellCt = '%d';
    size_Fr_ManualCellCt = [1 Inf]; % Inf can only be # Col

    Fr_ManualCellCt_Out = fscanf(fileID_Fr_ManualCellCt,formatSpec_Fr_ManualCellCt,size_Fr_ManualCellCt);
    Fr_ManualCellCt_Out = Fr_ManualCellCt_Out';

    fclose(fileID_Fr_ManualCellCt);

    if ZCropParam_In.Option == 'y'
        Fr_ManualCellCt_Out(Fr_ManualCellCt_Out>ZCropParam_In.ZEnd) = [];
        Fr_ManualCellCt_Out = Fr_ManualCellCt_Out-ZCropParam_In.ZStart+1;
        Fr_ManualCellCt_Out(Fr_ManualCellCt_Out<1) = [];
    end
else
    Fr_ManualCellCt_Out = [];
    warning('Manual CellCt information not available in directory.')
 end

clearvars fileID_Fr_ManualCellCt formatSpec_Fr_ManualCellCt size_Fr_ManualCellCt;
clearvars Fr_ManualCellCt_FilenameString;