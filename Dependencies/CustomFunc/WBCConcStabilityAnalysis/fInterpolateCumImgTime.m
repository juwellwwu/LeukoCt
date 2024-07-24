function [GroupID_CumVolvsCumTime_mCell_Out] = fInterpolateCumImgTime(GroupID_AllVidData_gCell_In,UI_In)

%% Designate interpolation step

CumVolIntp_Step_In = UI_In.CumVolIntp_Step;
fprintf('Interpolate Cumulative ImgTime (s) from Cumulative blood flow volume (nL) for each GroupID (Vessel identifier) ...\n');


%% Organize x (= Volume) ,y (= Time) data used for interpolation

x_gCell = cell(height(GroupID_AllVidData_gCell_In),1);
y_gCell = cell(height(GroupID_AllVidData_gCell_In),1);

for ugCt = 1:height(GroupID_AllVidData_gCell_In) 

    x_gCell{ugCt,1} = zeros(height(GroupID_AllVidData_gCell_In{ugCt,2}),3,2); 
    
    x_gCell{ugCt,1}(:,1,1) = GroupID_AllVidData_gCell_In{ugCt,2}(:,31); 
    x_gCell{ugCt,1}(:,1,2) = GroupID_AllVidData_gCell_In{ugCt,2}(:,34); 

    x_gCell{ugCt,1}(:,2,1) = GroupID_AllVidData_gCell_In{ugCt,2}(:,33); 
    x_gCell{ugCt,1}(:,2,2) = GroupID_AllVidData_gCell_In{ugCt,2}(:,36); 

    x_gCell{ugCt,1}(:,3,1) = GroupID_AllVidData_gCell_In{ugCt,2}(:,32); 
    x_gCell{ugCt,1}(:,3,2) = GroupID_AllVidData_gCell_In{ugCt,2}(:,35); 
    
    y_gCell{ugCt,1} = repmat(GroupID_AllVidData_gCell_In{ugCt,2}(:,3),1,1,2); 

end

clearvars ugCt;


%% Interpolate

GroupID_CumVolvsCumTime_mCell_Out = cell(2,1); 

for mCt = 1:2 

    CumVolMax_All = max(cellfun(@(x) max(x(:,1:3,:),[],'all'), x_gCell(~cellfun(@isempty,x_gCell))),[],1);
    CumVolMax_All = CumVolIntp_Step_In*ceil(CumVolMax_All/CumVolIntp_Step_In);

    % = Set up empty matrix that will hold interpolated data
    GroupID_CumVolvsCumTime = nan(height([0:CumVolIntp_Step_In:CumVolMax_All]'),1+height(x_gCell),3); 

    % = Load Col 1 (interpolated Cum imaging time (s) or blood volume (nL))
    GroupID_CumVolvsCumTime(:,1,:) = repmat([0:CumVolIntp_Step_In:CumVolMax_All]',1,1,3);

    % = Load Col 2-end (interpolated Cum WBC Conc (K/uL))
    for ugCt = 1:height(x_gCell) % 

        CumVolMax = max(x_gCell{ugCt,1}(:,1:3,mCt),[],'all');
        CumVolMax = CumVolIntp_Step_In*ceil(CumVolMax/CumVolIntp_Step_In);        
        CumVol_q = [0:CumVolIntp_Step_In:CumVolMax]'; 

        % i. k=1: Interpotation of Cum Img Time (y) by Cum Vol (x)       
        CumTime_q = interp1(vertcat(0,x_gCell{ugCt,1}(:,1,mCt)),vertcat(0,y_gCell{ugCt,1}(:,1,mCt)),CumVol_q); 
        GroupID_CumVolvsCumTime(1:height(CumTime_q),ugCt+1,1) = CumTime_q;

       % ii. k=2: Interpotation of Cum Img Time (y) by Cum Vol_Lo (x)
        CumTime_q = interp1(vertcat(0,x_gCell{ugCt,1}(:,2,mCt)),vertcat(0,y_gCell{ugCt,1}(:,1,mCt)),CumVol_q); 
        GroupID_CumVolvsCumTime(1:height(CumTime_q),ugCt+1,2) = CumTime_q;

        % iii. k=3: Interpotation of Cum Img Time (y) by Cum Vol_Hi (x)
        CumTime_q = interp1(vertcat(0,x_gCell{ugCt,1}(:,3,mCt)),vertcat(0,y_gCell{ugCt,1}(:,1,mCt)),CumVol_q); 
        GroupID_CumVolvsCumTime(1:height(CumTime_q),ugCt+1,3) = CumTime_q;

    end

    GroupID_CumVolvsCumTime_mCell_Out{mCt,1} = GroupID_CumVolvsCumTime;

end