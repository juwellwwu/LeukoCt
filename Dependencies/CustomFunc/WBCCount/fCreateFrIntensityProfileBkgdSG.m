function [LinePlot_Bkgd_Out] = fCreateFrIntensityProfileBkgdSG(LinePlot_In,sgolay_Order_In,sgolay_FrLength_In,Z_PxLength_In,PkRmvOption,PkRmvMethodOption) 

NumBlk_In = width(LinePlot_In);

%% Prepare Outlier (WBC removed) Intensity Profile for Background Preparation

if PkRmvOption == 'y'
    LinePlot_PeakRmv = fRmvOutlierPeak(LinePlot_In,Z_PxLength_In,PkRmvMethodOption);
else
    LinePlot_PeakRmv = LinePlot_In;
end

%% == Savitzky-Golay Smoothing

% Savitzky-Golay Smooth frame length for creating background 
% Must be odd
if rem(sgolay_FrLength_In,2)<eps
    sgolay_FrLength = sgolay_FrLength_In+1;
else
    sgolay_FrLength = sgolay_FrLength_In;
end

% = Smooth PeakRmv line with Savitzky-Golay filter
LinePlot_Bkgd_Out = zeros(size(LinePlot_In),'single');
for BlkCt = 1:NumBlk_In

    LinePlot_Bkgd_Out(:,BlkCt) = sgolayfilt(double(LinePlot_PeakRmv(:,BlkCt)),double(sgolay_Order_In),double(sgolay_FrLength));

    % Remove outliers at 2 ends, if present
    LinePlot_Bkgd_Out(1:ceil(1.00*sgolay_FrLength),BlkCt) = filloutliers(LinePlot_Bkgd_Out(1:ceil(1.00*sgolay_FrLength),BlkCt),'center',"median");
    LinePlot_Bkgd_Out((end-ceil(1.00*sgolay_FrLength)+1):end,BlkCt) = filloutliers(LinePlot_Bkgd_Out((end-ceil(1.00*sgolay_FrLength)+1):end,BlkCt),'center',"median");

end

