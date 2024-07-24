function [Skel_xy_Cell_Output,Skel_bdDist_Cell_Output,Skel_xy_Flag_Cell_Output] = ...
    fCreateSkelfromThreshPerim(ThreshPerim_xy_tpCell_Input,ThreshPerim_xy_Flag_tpCell_Input,Vessel_RadiusL_px_Input,NumSkel_VsWidth_Input, bdDist_Lim)

SkelBlockLength_Input = height(ThreshPerim_xy_tpCell_Input{1,1});
NumBlk_Input = size(ThreshPerim_xy_tpCell_Input,3);

% == Initiate cell for skeletons
Skel_xy_Cell_Output = repmat({zeros(SkelBlockLength_Input,2,'double')},NumSkel_VsWidth_Input,1,NumBlk_Input);

% == Initiate cell for recording distance between skeletons
Skel_bdDist_Cell_Output = repmat({zeros(SkelBlockLength_Input,2,'double')},1,1,NumBlk_Input); 

% == Initiate cell for flagging SkelCt, SkelPxCt too close to vessel boundary
Skel_xy_Flag_Cell_Output = repmat({zeros(SkelBlockLength_Input,1,'logical')},NumSkel_VsWidth_Input,1,NumBlk_Input);

% == Load original skeleton as 1st element
Skel_xy_Cell_Output(1,1,:) = ThreshPerim_xy_tpCell_Input(1,1,:); 

% == Create and load new skeletons

% = Find time segment of max radius for each SkelPxCt; Bd = boundary
[~,Vessel_RadiusL_TimeMaxIdx_Bd1] = max(Vessel_RadiusL_px_Input(:,1:2:end),[],2);
[~,Vessel_RadiusL_TimeMaxIdx_Bd2] = max(Vessel_RadiusL_px_Input(:,2:2:end),[],2);

% = Shift time segment index back to index of ThreshPerim_Crop_xy_Cell
Vessel_RadiusL_TimeMaxIdx_Bd1 = reshape(Vessel_RadiusL_TimeMaxIdx_Bd1*2,SkelBlockLength_Input,1,NumBlk_Input);
Vessel_RadiusL_TimeMaxIdx_Bd2 = reshape((Vessel_RadiusL_TimeMaxIdx_Bd2*2+1),SkelBlockLength_Input,1,NumBlk_Input);

% = Create new skeletons
for bCt = 1:2 
    for BlkCt = 1:NumBlk_Input

        if bCt == 1
            SkelIdx_Bd = (2:2:(NumSkel_VsWidth_Input-1)); 
            Vessel_RadiusL_TimeMaxIdx_Bd = Vessel_RadiusL_TimeMaxIdx_Bd1;
        elseif bCt ==2
            SkelIdx_Bd = (3:2:NumSkel_VsWidth_Input); 
            Vessel_RadiusL_TimeMaxIdx_Bd = Vessel_RadiusL_TimeMaxIdx_Bd2;
        end

        for SkelPxCt = 1:SkelBlockLength_Input

            % Load (x,y) of Img_Thresh Perimeter w/ max radius of all time points
            Skel_xy_Cell_Output{SkelIdx_Bd(end),1,BlkCt}(SkelPxCt,:) = ...
                ThreshPerim_xy_tpCell_Input{Vessel_RadiusL_TimeMaxIdx_Bd(SkelPxCt,1,BlkCt),1,BlkCt}(SkelPxCt,:);

            % Calculate new (x,y) to spread the new skeletons evenly between original skeleton and boundary
            xy_linspace = ... 
                horzcat((linspace(Skel_xy_Cell_Output{1,1,BlkCt}(SkelPxCt,1),Skel_xy_Cell_Output{SkelIdx_Bd(end),1,BlkCt}(SkelPxCt,1),0.5*(NumSkel_VsWidth_Input-1)+1))',...
                (linspace(Skel_xy_Cell_Output{1,1,BlkCt}(SkelPxCt,2),Skel_xy_Cell_Output{SkelIdx_Bd(end),1,BlkCt}(SkelPxCt,2),0.5*(NumSkel_VsWidth_Input-1)+1))');
            xy_linspace(1,:) = []; 

            % Load x,y into intermediate skeletons
            for SkelCt = 1:(numel(SkelIdx_Bd)-1)
                Skel_xy_Cell_Output{SkelIdx_Bd(SkelCt),1,BlkCt}(SkelPxCt,:) = round(xy_linspace(SkelCt,:));
            end
            
            % Load distance between skeletons into cell 
            Skel_bdDist_Cell_Output{1,BlkCt}(SkelPxCt,bCt) = pdist2(xy_linspace(2,:),xy_linspace(1,:),'euclidean');

            % Create flag cell  
            Skel_xy_Flag_Cell_Output{SkelIdx_Bd(SkelCt),1,BlkCt}(SkelPxCt,1) = ...
                ThreshPerim_xy_Flag_tpCell_Input{Vessel_RadiusL_TimeMaxIdx_Bd(SkelPxCt,1,BlkCt),1,BlkCt}(SkelPxCt,1);
  
            % i) Skel_xy too close to boundary (bdDist_Lim)
            xy_linspace_Dist = pdist2(xy_linspace(end,:),xy_linspace(:,:),'euclidean'); % distance between vessel boundary and (x,y)
            xy_linspace_Flag = zeros(numel(xy_linspace_Dist),1,'logical');
            xy_linspace_Flag(xy_linspace_Dist'<bdDist_Lim) = 1;

            for SkelCt = 1:numel(SkelIdx_Bd)
                Skel_xy_Flag_Cell_Output{SkelIdx_Bd(SkelCt),1,BlkCt}(SkelPxCt,1) = ...
                    xy_linspace_Flag(SkelCt,1) || ThreshPerim_xy_Flag_tpCell_Input{Vessel_RadiusL_TimeMaxIdx_Bd(SkelPxCt,1,BlkCt),1,BlkCt}(SkelPxCt,1);
            end

        end

    end
end