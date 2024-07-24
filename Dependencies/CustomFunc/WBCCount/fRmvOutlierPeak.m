function [LinePlot_OutlierRmv_Out] = fRmvOutlierPeak(LinePlot_In,Z_PxLength_In,BkgdOption) 

NumFr_1sec = 1/Z_PxLength_In;
NumBlk_In = width(LinePlot_In);

%% == Background Peak Removal

if BkgdOption == 1

    LinePlot_OutlierRmv_Out = zeros(size(LinePlot_In),'single');
    LinePlot_PeakRmv_UpLim = zeros(1,width(LinePlot_In),'single');
    LinePlot_PeakRmv_UpLimOld = zeros(1,width(LinePlot_In),'single');
    LinePlot_Bkgd_Out = zeros(size(LinePlot_In),'single');

    for BlkCt = 1:NumBlk_In

        LinePlot_OutlierRmv_Out(:,BlkCt) = LinePlot_In(:,BlkCt);

        LinePlot_PeakRmv_UpLim(1,BlkCt) = max(LinePlot_OutlierRmv_Out(:,BlkCt),[],'all'); % initiate
        LinePlot_PeakRmv_UpLimOld(1,BlkCt) = 9999; % initiate

        while_Ct=0;

        while (LinePlot_PeakRmv_UpLim(1,BlkCt)/LinePlot_PeakRmv_UpLimOld(1,BlkCt)<0.95)
            LinePlot_PeakRmv_UpLimOld(1,BlkCt) = LinePlot_PeakRmv_UpLim(1,BlkCt);
            LinePlot_PeakRmv_UpLim(1,BlkCt) = ...
                median(LinePlot_OutlierRmv_Out(:,BlkCt)) + 3.00*1.4826*mad(LinePlot_OutlierRmv_Out(:,BlkCt,1)); % pre-20230526: 3.00*1.4826
            LinePlot_OutlierRmv_Out(LinePlot_OutlierRmv_Out(:,BlkCt)>LinePlot_PeakRmv_UpLim(1,BlkCt),BlkCt) = ...
                median(LinePlot_OutlierRmv_Out(:,BlkCt));
            while_Ct = while_Ct+1;
        end

    end

end


%% == Background Peak Removal

if BkgdOption == 2

    LinePlot_movmin = zeros(size(LinePlot_In),'single');
    LinePlot_movmax= zeros(size(LinePlot_In),'single');
    LinePlot_movrng= zeros(size(LinePlot_In),'single');
    LinePlot_UpLim = zeros(size(LinePlot_In),'single');
    LinePlot_RoughBkgd = zeros(size(LinePlot_In),'single');

    for BlkCt = 1:NumBlk_In

        LinePlot_movmin(:,BlkCt) = movmin(LinePlot_In,NumFr_1sec);
        LinePlot_movmax(:,BlkCt) = movmax(LinePlot_In,NumFr_1sec);
        LinePlot_movrng(:,BlkCt) = LinePlot_movmax(:,BlkCt)-LinePlot_movmin(:,BlkCt);
        LinePlot_UpLim(:,BlkCt) = LinePlot_movmin(:,BlkCt) + prctile(LinePlot_movrng(:,BlkCt),10);

        LinePlot_RoughBkgd(:,BlkCt) = LinePlot_In(:,BlkCt);
        LinePlot_RoughBkgd(LinePlot_RoughBkgd(:,BlkCt)>LinePlot_UpLim(:,BlkCt),BlkCt) = ...
            LinePlot_UpLim(LinePlot_RoughBkgd(:,BlkCt)>LinePlot_UpLim(:,BlkCt),BlkCt);

    end

    % Check
    % figure;
    % plot(LinePlot_In(:,BlkCt));
    % hold on;
    % plot(LinePlot_RoughBkgd(:,BlkCt));
    
    % i) Create rough background w/ major peaks removed
    LinePlot_RoughPeakRmv = zeros(size(LinePlot_In),'single');
    for BlkCt = 1:NumBlk_In
        LinePlot_RoughBkgd(:,BlkCt) =  filloutliers(LinePlot_RoughBkgd(:,BlkCt),"linear",'median','ThresholdFactor',3.0);
        for smCt = (floor(height(LinePlot_In)/(NumFr_1sec*5)):-1:1) % sequential in 5 sec
            LinePlot_RoughBkgd(:,BlkCt) =  filloutliers(LinePlot_RoughBkgd(:,BlkCt),"linear",'movmedian',(NumFr_1sec*5)*smCt,'ThresholdFactor',3.0);
        end
        LinePlot_RoughPeakRmv(:,BlkCt) = LinePlot_In(:,BlkCt)-LinePlot_RoughBkgd(:,BlkCt);
    end

    % Check
    % figure;
    % plot(LinePlot_In(:,BlkCt));
    % hold on;
    % plot(LinePlot_RoughBkgd(:,BlkCt));
    % hold on;
    % plot(LinePlot_RoughPeakRmv(:,BlkCt));

    % ii) Determine location of peaks on LinePlot_RoughPeakRmv
    LinePlot_NaN = zeros(size(LinePlot_In),'logical');

    for BlkCt = 1:NumBlk_In

        [pk_height,pk_loc,pk_width,pk_prom,pk_width_ext] = findpeaks_MOD(LinePlot_RoughPeakRmv,'WidthReference','halfprom');

        pk_loc(pk_height<0.1*max(pk_height,[],'all')) = [];
        pk_width(pk_height<0.1*max(pk_height,[],'all')) = [];
        pk_prom(pk_height<0.1*max(pk_height,[],'all')) = [];
        pk_width_ext(pk_height<0.1*max(pk_height,[],'all'),:) = [];
        pk_height(pk_height<0.1*max(pk_height,[],'all')) = []; 

        for pkCt = 1:numel(pk_height)

            pk_StartIdx = pk_loc(pkCt) - round(1.7*(pk_loc(pkCt)-pk_width_ext(pkCt,1)));
            pk_StartIdx(pk_StartIdx<1) = 1;
            pk_EndIdx = pk_loc(pkCt) + round(1.7*(pk_width_ext(pkCt,2)-pk_loc(pkCt)));
            pk_EndIdx(pk_EndIdx>height(LinePlot_In)) = height(LinePlot_In);

            LinePlot_NaN(pk_StartIdx:pk_EndIdx,BlkCt) = 1;

        end

    end

    % iii) Create refined peak-removed LinePlot_In for filtering
    LinePlot_OutlierRmv_Out = LinePlot_In;

    for BlkCt = 1:NumBlk_In
        LinePlot_OutlierRmv_Out(LinePlot_NaN(:,BlkCt),1) = NaN;
        LinePlot_OutlierRmv_Out(:,BlkCt) = fillmissing(LinePlot_OutlierRmv_Out(:,BlkCt),'linear');
    end

end