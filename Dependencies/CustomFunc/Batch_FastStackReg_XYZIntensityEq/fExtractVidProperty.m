function [SampleIDString_Out, XY_PxLength_Out, Z_PxLength_Out,PreReg_XY_PxLength_Out, Reg_XY_PxLength_Out, PostReg_XY_PxLength_Out] = ...
    fExtractVidProperty(VidProp_FilenameString_In,Vid_FilenameString_In, expr_SampleID_In)

% Extract SampleID from video name
[StartIdx,EndIdx] = regexp(Vid_FilenameString_In,expr_SampleID_In);
SampleIDString_Out = Vid_FilenameString_In(StartIdx:EndIdx);

% Remove Slice Information when search for video properties
[StartIdx,EndIdx] = regexp(SampleIDString_Out,'_S\d+-\d+');
if ~isempty(StartIdx)
    SampleIDString_NoSlice = SampleIDString_Out(1:(StartIdx-1));
else
    SampleIDString_NoSlice = SampleIDString_Out;
end

VidProp = readtable(VidProp_FilenameString_In); 

% Find Row idx matching Sample ID in Excel worksheet
SampleID_RowIdx = ... 
    find(cellfun(@(x) isequal(SampleIDString_NoSlice,x), VidProp.SampleID, 'UniformOutput', true)); 

% Extract pixel size information
XY_PxLength_Out = 1/VidProp.Recp_XYPxLength(SampleID_RowIdx); 
Z_PxLength_Out = 1/VidProp.Recp_ZPxLength(SampleID_RowIdx); 

PreReg_XY_PxLength_Out  = 1/VidProp.Recp_PreRegXYPxLength(SampleID_RowIdx);
Reg_XY_PxLength_Out = 1/VidProp.Recp_RegXYPxLength(SampleID_RowIdx);
PostReg_XY_PxLength_Out = 1/VidProp.Recp_PostRegXYPxLength(SampleID_RowIdx);

clearvars StartIdx EndIdx SampleID_RowIdx SampleIDString_NoSlice;