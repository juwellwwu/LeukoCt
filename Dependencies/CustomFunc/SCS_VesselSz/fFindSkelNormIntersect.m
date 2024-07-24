function [Intersect_DilateSkelPxIdx_SkCell_Output,Intersect_DilateSkel_xy_SkCell_Output,Img_Norm_Output,rp_Norm_Output] = ...
    fFindSkelNormIntersect(Img_Norm_Input,Skel_xy_Ord_Cell_Input,SkelPxCt_Input,Intersect_DilateSkelPxIdx_SkCell_Input,Intersect_DilateSkel_xy_SkCell_Input)

% = Keep only the chunk of normal line that contains the original
% skeleton pixel for which the normal is drawn
rp_Norm = regionprops(Img_Norm_Input,'Area','PixelList','PixelIdxList');

ismember_FlagID = zeros(height(rp_Norm),1,'logical'); 
for nsCt = 1:height(rp_Norm) 
    ismember_FlagID(nsCt,1) = ~ismember(Skel_xy_Ord_Cell_Input{1,1}(SkelPxCt_Input,:),rp_Norm(nsCt).PixelList,'rows');
end

if height(rp_Norm) == sum(ismember_FlagID)
    for nsCt = 1:height(rp_Norm) 
        [pDist,~] = pdist2(rp_Norm(nsCt).PixelList,Skel_xy_Ord_Cell_Input{1,1}(SkelPxCt_Input,:),@fMeasureChessboardDist,'Smallest',2);
        pDist(pDist>1) = NaN;
        if all(isnan(pDist))
            ismember_FlagID(nsCt,1) = 1;
        else
            ismember_FlagID(nsCt,1) = 0;
        end
    end
end

% = Update Normal image and regionprops
Img_Norm_Output = Img_Norm_Input;
Img_Norm_Output(cat(1,rp_Norm(ismember_FlagID).PixelIdxList)) = 0;

rp_Norm(ismember_FlagID) = [];
rp_Norm_Output = rp_Norm;

% = OPTIONAL: Check
% Build Skeleton Img
% % % SkelIdx = sub2ind(size(Img_Norm_Input),Skel_xy_Ord_Cell_Input{1,1}(:,2),Skel_xy_Ord_Cell_Input{1,1}(:,1));
% % % Img_Skel = zeros(size(Img_Norm_Input),'logical');
% % % Img_Skel(SkelIdx) = 1;
% % % figure; imshowpair(Img_Skel,Img_Norm_Output);

% = Determine where original skeleton's line normal intersects all
% skeletons

Intersect_DilateSkel_xy_SkCell_Output = Intersect_DilateSkel_xy_SkCell_Input;
Intersect_DilateSkelPxIdx_SkCell_Output = Intersect_DilateSkelPxIdx_SkCell_Input;

if ~isempty(rp_Norm)
    for SkelCt = 1:height(Skel_xy_Ord_Cell_Input)
        
        % = Simplest case: Skeleton normal intersects 
        [Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input},Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input},~] = ...
            intersect([Skel_xy_Ord_Cell_Input{SkelCt,1}(:,1) Skel_xy_Ord_Cell_Input{SkelCt,1}(:,2)], [rp_Norm.PixelList(:,1) rp_Norm.PixelList(:,2)],'rows');
        
        % = Next case: Skeleton normal do not intersect due to diagnonal line crossing but pixels not overlapping
        if isempty(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input})
            Img_Norm_Temp = imdilate(Img_Norm_Output,strel('disk',1));
            rp_Norm_Temp = regionprops(Img_Norm_Temp,'Area','PixelList','PixelIdxList');
            [Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input},Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input},~] = ...
                intersect([Skel_xy_Ord_Cell_Input{SkelCt,1}(:,1) Skel_xy_Ord_Cell_Input{SkelCt,1}(:,2)], [rp_Norm_Temp.PixelList(:,1) rp_Norm_Temp.PixelList(:,2)],'rows');
        end

        % = Clean up: if more than 1 intersection from simpliest or next case
        if (numel(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input})>1) 
            [Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input},SortIdx] = sort(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input},'ascend');
            Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input} = Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input}(SortIdx,:);
            if (max(abs(diff(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input})))==1) 
                Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input} = ... 
                    Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input}(ceil(0.5*height(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input})),:);
                Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input} = ... 
                    Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input}(ceil(0.5*height(Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input})),:);
            else
                [~,MinIdx] = ...
                    min(abs(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input}-repelem(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input-1},numel(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input}),1)));
                Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input} = ...
                    Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input}(MinIdx,:);
                Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input} = ... 
                    Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input}(MinIdx,:);
            end
        end

        if isempty(Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input})
            Intersect_DilateSkelPxIdx_SkCell_Output{SkelCt,SkelPxCt_Input} = NaN;
            Intersect_DilateSkel_xy_SkCell_Output{SkelCt,SkelPxCt_Input} = [-9999 -9999];
        end

    end

end

end
