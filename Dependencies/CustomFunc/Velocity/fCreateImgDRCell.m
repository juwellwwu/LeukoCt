function [Img_DR_Cell_Output] = ...
    fCreateImgDRCell(ImgStack_Input,Skel_xy_Cell_Input,NumSkel_Input,NumSeg_Input,NumBlk_Input,TimeParam_Input)

fprintf('Performing Digital Line Scan: Create cell of space-time diagrams...\n');

SkelBlockLength_Input = height(Skel_xy_Cell_Input{1,1});

Fr_Start_Input = TimeParam_Input.tSeg_FrStart;
Fr_End_Input = TimeParam_Input.tSeg_FrEnd;

Seg_FrLength_Input = TimeParam_Input.tSeg_FrLength;

Img_DR_Cell_Output = cell(NumSkel_Input,NumSeg_Input,NumBlk_Input);

for SegCt = 1:NumSeg_Input
    for SkelCt = 1:NumSkel_Input 
        for BlkCt = 1:NumBlk_Input

            Img_DR_Cell_Output{SkelCt,SegCt,BlkCt} = zeros(Seg_FrLength_Input,SkelBlockLength_Input,'single');
            for SkelPxCt=1:size(Skel_xy_Cell_Input{SkelCt,1,BlkCt},1) 
                Img_DR_Cell_Output{SkelCt,SegCt,BlkCt}(:,SkelPxCt) = ...
                    reshape(ImgStack_Input(Skel_xy_Cell_Input{SkelCt,1,BlkCt}(SkelPxCt,2),Skel_xy_Cell_Input{SkelCt,1,BlkCt}(SkelPxCt,1),Fr_Start_Input(SegCt):Fr_End_Input(SegCt)),Seg_FrLength_Input,1);
            end

        end
    end
end