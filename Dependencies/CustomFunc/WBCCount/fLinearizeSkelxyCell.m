function [Skel_LinIdx_Mtx_Out,Skel_x_Mtx_Out,Skel_y_Mtx_Out] = fLinearizeSkelxyCell(ImgStack_In,Skel_xy_Cell_In)

NumBlk_In = size(Skel_xy_Cell_In,3);

sub2ind_wrapper = @(x) sub2ind([height(ImgStack_In) width(ImgStack_In)],x(:,2),x(:,1));

% Create Skel_LinIdx_Mtx, in which
% # Row: SkelPxCt for all vessel blocks (SkelBlockLength*NumBlk total)
% # Col: Skel Idx (NumSkel_VsWidth total)
Skel_LinIdx_Mtx_Out = cellfun(sub2ind_wrapper,Skel_xy_Cell_In,'UniformOutput',false);
Skel_LinIdx_Mtx_Out = cell2mat(permute(Skel_LinIdx_Mtx_Out,[2 1 3])); 
Skel_LinIdx_Mtx_Out = mat2cell(Skel_LinIdx_Mtx_Out,height(Skel_LinIdx_Mtx_Out),width(Skel_LinIdx_Mtx_Out),repelem(1,NumBlk_In));
Skel_LinIdx_Mtx_Out = cell2mat(permute(Skel_LinIdx_Mtx_Out,[3 1 2]));

% In similar organization, create Skel_x_Mtx, and Skel_y_Mtx
[Skel_x_Mtx_Out,Skel_y_Mtx_Out] = ind2sub([height(ImgStack_In) width(ImgStack_In)],Skel_LinIdx_Mtx_Out);
