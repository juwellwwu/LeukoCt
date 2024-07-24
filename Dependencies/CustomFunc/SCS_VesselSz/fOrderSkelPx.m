function SkelPxList_Ord = fOrderSkelPx(SkelPxList,Endpt_Loc)

[pDist,pDist_Loc] = pdist2(SkelPxList,SkelPxList,@fMeasureChessboardDist,'Smallest',height(SkelPxList)); 
pDist = pDist';
pDist = pDist(:,2:end);
pDist_Loc = pDist_Loc';
pDist_Loc = pDist_Loc(:,2:end);
pDist_Loc(pDist>1) = NaN;
pDist(:,find(all(isnan(pDist_Loc),1))) = []; 
pDist_Loc(:,find(all(isnan(pDist_Loc),1))) = []; 

Skel_Order = zeros(size(SkelPxList,1),1,'single'); 
Skel_Order_Loc = zeros(size(SkelPxList,1),1,'single');

Skel_Order_Loc(1,1) = Endpt_Loc; 

for Skel_Order_PxCt = 1:(size(SkelPxList,1)-1)

    if Skel_Order_PxCt==1

        Skel_Order(Skel_Order_Loc(Skel_Order_PxCt,1),1) = Skel_Order_PxCt;
        Skel_Order_Loc(Skel_Order_PxCt+1,1) =  pDist_Loc(Skel_Order_Loc(Skel_Order_PxCt,1),1);

    elseif Skel_Order_PxCt>1  

        if find(~isnan(pDist_Loc(Skel_Order_Loc(Skel_Order_PxCt,1),:))) > 0
            Skel_Order_Loc(Skel_Order_PxCt+1,1) =  pDist_Loc(Skel_Order_Loc(Skel_Order_PxCt,1),min(find(~isnan(pDist_Loc(Skel_Order_Loc(Skel_Order_PxCt,1),:)))));
        end

    end

    pDist_Loc(pDist_Loc==Skel_Order_Loc(Skel_Order_PxCt,1)) = NaN; 

end

% Reorder Skeleton pixels in order specified Skel_Order_Org
SkelPxList_Ord = SkelPxList(Skel_Order_Loc,:);

