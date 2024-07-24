function [WinSkelPxIdx_Cell_Out,WinArea_px_Out] = fCalcSkelPxIdxGateXYWinArea(ImgStack_In,SCS_Param_In,WBCtoWinAreaRatio_In,XY_PxLength_In)

Skel_xy_Cell_In = SCS_Param_In.Skel_xy_Cell;

%% Window area 

WBCArea_px = pi*(0.5*12.0/XY_PxLength_In).^2;
WinArea_px_Out = WBCArea_px.*1./WBCtoWinAreaRatio_In;


%% Re-organize Skel (x,y) 

[Skel_LinIdx_Mtx,~,~] = fLinearizeSkelxyCell(ImgStack_In,Skel_xy_Cell_In);


%% Find area enclosed by pixel fill of each SkelPxCt "step"

SkelPxFillMaskArea = zeros(height(Skel_LinIdx_Mtx),1,'single');
SkelPxFillMask = zeros(height(ImgStack_In),width(ImgStack_In),height(Skel_LinIdx_Mtx),'logical');

NumSkel_PxWinAreaCalc = floor(0.5*0.8*width(Skel_LinIdx_Mtx))*2+1; 

for SkelPxCt = 1:height(Skel_LinIdx_Mtx)
    if SkelPxCt>1
        [~,SkelPxFillMaskArea(SkelPxCt,1)] = fConvertSkelxytoMask(ImgStack_In,Skel_LinIdx_Mtx((SkelPxCt-1):SkelPxCt,1:NumSkel_PxWinAreaCalc),'y');
    end
end


%% OPTIONAL: Check that PxFill of each SkelPxCt "step" adds back to Vessel ROI

% % % cmap = colormap(turbo(height(Skel_LinIdx_Mtx)));
% % % figure;
% % % imshow(ind2rgb(max(uint16(SkelPxFillMask).*uint16(permute((1:1:height(Skel_LinIdx_Mtx)),[1 3 2])),[],3),cmap),'InitialMagnification',200);
% % % 
% % % clearvars cmap;


%% Find closest Window Length, in SkelPx, that corresponds to Win_PxArea_Init

WinSkelPxIdx = zeros(height(Skel_LinIdx_Mtx),numel(WinArea_px_Out));

for SkelPxCt = 1:height(Skel_LinIdx_Mtx)
     
    SkelPxFillCumArea = SkelPxFillMaskArea;
    SkelPxFillCumArea(1:SkelPxCt,1) = 0;

    SkelPxFillCumArea = cumsum(SkelPxFillCumArea,1);
    [pDist, pDist_Idx] = pdist2(SkelPxFillCumArea,single(WinArea_px_Out'),'euclidean','Smallest',4); 
    
    for aCt = 1:numel(WinArea_px_Out)
        if pDist(1,aCt)>mean(SkelPxFillMaskArea(pDist_Idx(2:end,aCt),1))
            pDist_Idx(1,aCt) = NaN;
        end
    end
    
    WinSkelPxIdx(SkelPxCt,:,2) = pDist_Idx(1,:);
    WinSkelPxIdx(SkelPxCt,:,1) = SkelPxCt;

end

WinSkelPxIdx_Cell_Out = mat2cell(WinSkelPxIdx,height(Skel_LinIdx_Mtx),repelem(1,numel(WinArea_px_Out)),2);
WinSkelPxIdx_Cell_Out = cellfun(@(x) permute(x,[3 1 2])',WinSkelPxIdx_Cell_Out,'UniformOutput',false);

clearvars SkelPx_WinPxArea_StartEndIdx;


%%

close all;

clearvars pDist pDist_Idx SkelPxFillCumArea;
clearvars Skel_LinIdx_Mtx;
clearvars SkelPxFillArea;
clearvars aCt SkelPxCt;
