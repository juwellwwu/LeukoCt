function [Skel_xy_Ord_Cell_Output,Skel_xy_Flag_Cell_Output] = ...
    fDefineVesselBlock(NumBlk_Input,SkelBlockLength_Input,Skel_xy_Ord_Cell_Input,Skel_xy_Flag_Cell_Input)

if NumBlk_Input<1
    error('Selected Vessel ROI is of insufficient length for analysis. Terminate.');
end

Skel_xy_Ord_Cell_Output = Skel_xy_Ord_Cell_Input;
Skel_xy_Flag_Cell_Output = Skel_xy_Flag_Cell_Input;

for SkelCt = 1:height(Skel_xy_Ord_Cell_Input)
    if size(Skel_xy_Ord_Cell_Output{SkelCt,1},1) > (SkelBlockLength_Input*NumBlk_Input+1) % Trim off remainder pixel at both ends of skeleton

        dim1Dist = horzcat(ceil(0.5*(size(Skel_xy_Ord_Cell_Output{SkelCt,1},1)-SkelBlockLength_Input*NumBlk_Input)),...
            repmat(SkelBlockLength_Input,1,(NumBlk_Input-1)),size(Skel_xy_Ord_Cell_Output{SkelCt,1},1)-(size(Skel_xy_Ord_Cell_Output{SkelCt,1},1)-SkelBlockLength_Input*NumBlk_Input)-SkelBlockLength_Input*(NumBlk_Input-1),...
            floor(0.5*(size(Skel_xy_Ord_Cell_Output{SkelCt,1},1)-SkelBlockLength_Input*NumBlk_Input)));

        Skel_xy_Ord_Cell_Output{SkelCt,1}  = mat2cell(Skel_xy_Ord_Cell_Output{SkelCt,1},dim1Dist,2);
        Skel_xy_Ord_Cell_Output{SkelCt,1}(1) = [];
        Skel_xy_Ord_Cell_Output{SkelCt,1}(end) = [];
        
        Skel_xy_Flag_Cell_Output{SkelCt,1}  = mat2cell(Skel_xy_Flag_Cell_Output{SkelCt,1},dim1Dist,1);
        Skel_xy_Flag_Cell_Output{SkelCt,1}(1) = [];
        Skel_xy_Flag_Cell_Output{SkelCt,1}(end) = [];

    elseif (size(Skel_xy_Ord_Cell_Output{SkelCt,1},1) - (SkelBlockLength_Input*NumBlk_Input)) == 1 % Trim off one pixel at end of skeleton

        dim1Dist = horzcat(repmat(SkelBlockLength_Input,1,(NumBlk_Input-1)),size(Skel_xy_Ord_Cell_Output{SkelCt,1},1)-1-SkelBlockLength_Input*(NumBlk_Input-1),1);

        Skel_xy_Ord_Cell_Output{SkelCt,1}  = mat2cell(Skel_xy_Ord_Cell_Output{SkelCt,1},dim1Dist,2);
        Skel_xy_Ord_Cell_Output{SkelCt,1}(end) = [];

        Skel_xy_Flag_Cell_Output{SkelCt,1}  = mat2cell(Skel_xy_Flag_Cell_Output{SkelCt,1},dim1Dist,1);
        Skel_xy_Flag_Cell_Output{SkelCt,1}(end) = [];

    else  

        dim1Dist = horzcat(repmat(SkelBlockLength_Input,1,(NumBlk_Input-1)), size(Skel_xy_Ord_Cell_Output{SkelCt,1},1)-SkelBlockLength_Input*(NumBlk_Input-1));

        Skel_xy_Ord_Cell_Output{SkelCt,1}  = mat2cell(Skel_xy_Ord_Cell_Output{SkelCt,1},dim1Dist,2);
        Skel_xy_Flag_Cell_Output{SkelCt,1}  = mat2cell(Skel_xy_Flag_Cell_Output{SkelCt,1},dim1Dist,1);

    end
end

Skel_xy_Ord_Cell_Output = repmat(Skel_xy_Ord_Cell_Output,1,1,NumBlk_Input);
Skel_xy_Flag_Cell_Output = repmat(Skel_xy_Flag_Cell_Output,1,1,NumBlk_Input);

for SkelCt = 1:height(Skel_xy_Ord_Cell_Input)
    for BlkCt = 1:NumBlk_Input
        Skel_xy_Ord_Cell_Output{SkelCt,1,BlkCt} = cell2mat(Skel_xy_Ord_Cell_Output{SkelCt,1,BlkCt}(BlkCt,1));
        Skel_xy_Flag_Cell_Output{SkelCt,1,BlkCt} = cell2mat(Skel_xy_Flag_Cell_Output{SkelCt,1,BlkCt}(BlkCt,1));
    end
end