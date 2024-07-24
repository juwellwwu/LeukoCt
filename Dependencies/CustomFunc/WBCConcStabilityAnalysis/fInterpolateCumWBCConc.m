function [GroupID_CumTimeorVolvsCumConc_mCell_Out] = fInterpolateCumWBCConc(GroupID_AllVidData_gCell_In,TimeorVolOption_In,UI_In)

%% Elect interpolation step based on whether to interpolate against cumulative imaging time or volume.

if TimeorVolOption_In == 1 
    CumTorVIntp_Step_In = UI_In.CumTimeIntp_Step; 
    fprintf('Interpolate Cumulative WBC Concentration (K/uL) vs (cumulative) imaging time (s) for each GroupID (Vessl identifier) ...\n');
elseif TimeorVolOption_In == 2
    CumTorVIntp_Step_In = UI_In.CumVolIntp_Step;
    fprintf('Interpolate Cumulative WBC Concentration (K/uL) vs cumulative blood flow volume (nL) for each GroupID (Vessel identifier) ...\n');
end


%% Organize x (= Time or Volume) ,y (= WBC Conc) data used for interpolation

x_gCell = cell(height(GroupID_AllVidData_gCell_In),1);
y_gCell = cell(height(GroupID_AllVidData_gCell_In),1);

for ugCt = 1:height(GroupID_AllVidData_gCell_In)
    
    x_gCell{ugCt,1} = zeros(height(GroupID_AllVidData_gCell_In{ugCt,2}),1,2); 
    
    if TimeorVolOption_In == 1 
        x_gCell{ugCt,1} = repmat(GroupID_AllVidData_gCell_In{ugCt,2}(:,3),1,1,2); 
    elseif TimeorVolOption_In == 2 
        x_gCell{ugCt,1}(:,:,1) = GroupID_AllVidData_gCell_In{ugCt,2}(:,31); 
        x_gCell{ugCt,1}(:,:,2) = GroupID_AllVidData_gCell_In{ugCt,2}(:,34); 
    end

    y_gCell{ugCt,1} = zeros(height(GroupID_AllVidData_gCell_In{ugCt,2}),3,2); 
    y_gCell{ugCt,1}(:,1,1) = GroupID_AllVidData_gCell_In{ugCt,2}(:,37); 
    y_gCell{ugCt,1}(:,1,2) = GroupID_AllVidData_gCell_In{ugCt,2}(:,40); 

    y_gCell{ugCt,1}(:,2,1) = GroupID_AllVidData_gCell_In{ugCt,2}(:,38); 
    y_gCell{ugCt,1}(:,2,2) = GroupID_AllVidData_gCell_In{ugCt,2}(:,41); 

    y_gCell{ugCt,1}(:,3,1) = GroupID_AllVidData_gCell_In{ugCt,2}(:,39);
    y_gCell{ugCt,1}(:,3,2) = GroupID_AllVidData_gCell_In{ugCt,2}(:,42);

end

clearvars ugCt;


%% Interpolate

GroupID_CumTimeorVolvsCumConc_mCell_Out = cell(2,1); 

for mCt = 1:2 

    CumTorVMax_All = max(cellfun(@(x) max(x(:,1,:),[],'all'), x_gCell(~cellfun(@isempty,x_gCell))),[],1);

    % = Set up empty matrix that will hold interpolated data
    GroupID_CumTorVvsCumConc = nan(floor(CumTorVMax_All/CumTorVIntp_Step_In)+1,1+height(x_gCell),3);

    % = Load Col 1 (interpolated Cum imaging time (s) or blood volume (nL))
    GroupID_CumTorVvsCumConc(:,1,:) = repmat([0:CumTorVIntp_Step_In:CumTorVMax_All]',1,1,3);

    % = Load Col 2-end (interpolated Cum WBC Conc (K/uL))
    for ugCt = 1:height(x_gCell) % # Group ID

        if ~isempty(x_gCell{ugCt,1}) 

            CumTorVMax = max(x_gCell{ugCt,1}(:,1,:),[],'all');
            CumTorV_q = [0:CumTorVIntp_Step_In:(CumTorVIntp_Step_In*floor(CumTorVMax/CumTorVIntp_Step_In))]'; 

            % k=1; Cum WBC Conc
            CumConc_q = interp1(vertcat(0,x_gCell{ugCt,1}(:,1,mCt)),vertcat(0,y_gCell{ugCt,1}(:,1,mCt)),CumTorV_q); 
            GroupID_CumTorVvsCumConc(1:height(CumConc_q),ugCt+1,1) = CumConc_q;

            % k=2; Cum WBC_Conc_Lo based on Cum Vol Mean+StdEr
            CumConc_q = interp1(vertcat(0,x_gCell{ugCt,1}(:,1,mCt)),vertcat(0,y_gCell{ugCt,1}(:,2,mCt)),CumTorV_q); 
            GroupID_CumTorVvsCumConc(1:height(CumConc_q),ugCt+1,2) = CumConc_q;

            % k=3; Cum WBC Conc based on Cum Vol Mean-StdEr
            CumConc_q = interp1(vertcat(0,x_gCell{ugCt,1}(:,1,mCt)),vertcat(0,y_gCell{ugCt,1}(:,3,mCt)),CumTorV_q); 
            GroupID_CumTorVvsCumConc(1:height(CumConc_q),ugCt+1,3) = CumConc_q;

        end

    end

    GroupID_CumTimeorVolvsCumConc_mCell_Out{mCt,1} = GroupID_CumTorVvsCumConc;

end
