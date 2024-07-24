function [Vessel_RadiusL_px_Output,Vessel_Diameter_px_Output] = ...
    fCalcVesselSize(ThreshPerim_xy_OrdBlk_Cell_Input)

NumBlk_Input = size(ThreshPerim_xy_OrdBlk_Cell_Input,3);
ThreshPerimBlockLength_Input = height(ThreshPerim_xy_OrdBlk_Cell_Input{1,1}); % = SkelPxCt

% = Calculate Vessel Radius
Vessel_RadiusL_px_Output = zeros(ThreshPerimBlockLength_Input,height(ThreshPerim_xy_OrdBlk_Cell_Input),NumBlk_Input);

for BlkCt = 1:NumBlk_Input
    for SkelPxCt = 1:ThreshPerimBlockLength_Input

        TP_xy_Temp = cell2mat(ThreshPerim_xy_OrdBlk_Cell_Input(:,1,BlkCt));
        TP_xy_Temp = TP_xy_Temp(SkelPxCt:ThreshPerimBlockLength_Input:end,:);
        
        % Distance between skeleton pixel and thresh perim
        Vessel_RadiusL_px_Output(SkelPxCt,:,BlkCt) = ...
            pdist2(ThreshPerim_xy_OrdBlk_Cell_Input{1,1,BlkCt}(SkelPxCt,:),TP_xy_Temp,'euclidean');
    end
end

Vessel_RadiusL_px_Output = ... 
    mat2cell(Vessel_RadiusL_px_Output,height(Vessel_RadiusL_px_Output),...
    width(Vessel_RadiusL_px_Output),repelem(1,NumBlk_Input));
Vessel_RadiusL_px_Output = cell2mat(permute(Vessel_RadiusL_px_Output,[3,1,2]));
Vessel_RadiusL_px_Output(:,1) = [];

% = Calculate Vessel Diameter
Vessel_Diameter_px_Output = ...
    Vessel_RadiusL_px_Output(:,1:2:end)+Vessel_RadiusL_px_Output(:,2:2:end);

% = Set Output Radius
Vessel_RadiusL_px_Output = max(Vessel_RadiusL_px_Output(:,1:2:end),Vessel_RadiusL_px_Output(:,2:2:end));
