function [Img_FFT_Filtered_Clean_Cell_Out,rp_FFT_Filtered_Clean_Cell_Out] = fCleanImgFFTFilteredCell(Img_FFT_Filtered_Cell_In,Img_Postprocess_Cell_In,XY_PxLength_In)

NumSkel_In = size(Img_FFT_Filtered_Cell_In,1);
NumSeg_In = size(Img_FFT_Filtered_Cell_In,2);
NumBlk_In = size(Img_FFT_Filtered_Cell_In,3);

SkelBlockLength_In = size(Img_FFT_Filtered_Cell_In{1,1,1},2);

rp_FFT_Filtered_Clean_Cell_Out = cell(NumSkel_In,NumSeg_In,NumBlk_In);
rp_FFT_Filtered_Clean_FlagBlobIdx_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In);
rp_FFT_Filtered_Clean_FlagBlobIdx2_Cell = cell(NumSkel_In,NumSeg_In,NumBlk_In); % 20240217: check if frame rate sufficient for velocity measurement

% Set minimum slope line length
SlopeLineLength_LoLimit = 30/XY_PxLength_In;
if SlopeLineLength_LoLimit<30/XY_PxLength_In
    SlopeLineLength_LoLimit = SkelBlockLength_In*0.75;
    warning(strcat('Due to short SkelBlockLength, SlopeLineLength_LoLimit is < 30um =',32, num2str(SlopeLineLength_LoLimit.*XY_PxLength_In,'%0.2f'), 'um.'));
end

MaxFeretMaxMinRatio_Cell = repmat({NaN},NumSkel_In,NumSeg_In,NumBlk_In); % cell w/ all elements = NaN

tic
for SegCt = 1:NumSeg_In 
    for SkelCt = 1:NumSkel_In
        for BlkCt = 1:NumBlk_In

            % = Clear residual bkgd
            Img_FFT_Filtered_Cell_In{SkelCt,SegCt,BlkCt}(Img_FFT_Filtered_Cell_In{SkelCt,SegCt,BlkCt}<0.015) = 0;

            % = Binarize
            Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = Img_FFT_Filtered_Cell_In{SkelCt,SegCt,BlkCt};

            if max(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},[],'all')>eps

                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = logical(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt});
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = bwmorph(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},'hbreak');
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = bwskel(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},'MinBranchLength',floor(SlopeLineLength_LoLimit));

                % = Quick Clean up (Feret Diameters time-costly to calculate)
                % 1: remove v. small blobs
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = ...
                    regionprops(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},'Area','PixelIdxList');
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.Area)<(0.80*SlopeLineLength_LoLimit)).PixelIdxList)) = 0;
                if max(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},[],'all')<eps % All blobs cleared
                    rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = [];
                    Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = single(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt});
                    continue;
                end

                % 2: remove + orientation (wrong flow direction)
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = ... 
                    regionprops(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},'Orientation','PixelIdxList');
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.Orientation)>=(0.00)).PixelIdxList)) = 0;
                if max(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},[],'all')<eps % All blobs cleared
                    rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = [];
                    Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = single(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt});
                    continue;
                end

                % 3: remove orientation w/in (3/SkelBlockLength)*180/pi deg of 0 or 5 deg of 90
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = ... % Update
                    regionprops(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},'Orientation','PixelIdxList','PixelList');
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(abs(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.Orientation))<(3/SkelBlockLength_In*180/pi)).PixelIdxList)) = 0;
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(abs(abs(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.Orientation))-90)<(5.00)).PixelIdxList)) = 0;
                if max(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},[],'all')<eps % All blobs cleared
                    rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = [];
                    Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = single(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt});
                    continue;
                end

                % Update regionprops after Quick Clean Up
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = ...
                    regionprops(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},'Area','PixelIdxList','PixelList','MaxFeretProperties','MinFeretProperties','Orientation');

            else

                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = zeros(size(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}),'logical');
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = [];
                MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt} = NaN;

            end

            % = Thorough Clean Up
            if ~isempty(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})

                rp_FFT_Filtered_Clean_FlagBlobIdx_Cell{SkelCt,SegCt,BlkCt} = ...
                    (cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MinFeretDiameter)>2.0+floor(0.05*cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MaxFeretDiameter))) | ...
                    (cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.Area)./cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MaxFeretDiameter)>1.10) | ...
                    (cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MaxFeretDiameter)<SlopeLineLength_LoLimit) | ...
                    (cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MaxFeretDiameter)<0.33*max(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MaxFeretDiameter)));

                % Set Flagged blobs intensity to 0
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(rp_FFT_Filtered_Clean_FlagBlobIdx_Cell{SkelCt,SegCt,BlkCt}).PixelIdxList)) = 0;
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = single(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt});

                % Remove Flagged blobs from regionprops
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(rp_FFT_Filtered_Clean_FlagBlobIdx_Cell{SkelCt,SegCt,BlkCt}) = [];

            end

            % = Ensure imaging frame rate sufficient for flow velocity
            if ~isempty(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})

                rp_y_range = zeros(height(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}),1,'single');
                for bCt = 1:height(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})
                    rp_y_range_Temp = rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(bCt).PixelList;
                    rp_y_range_Temp = range(rp_y_range_Temp(:,2));
                    rp_y_range(bCt,1) = rp_y_range_Temp;
                end

                rp_FFT_Filtered_Clean_FlagBlobIdx2_Cell{SkelCt,SegCt,BlkCt} = rp_y_range < 3;

                % Set Flagged blobs intensity to 0
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(rp_FFT_Filtered_Clean_FlagBlobIdx2_Cell{SkelCt,SegCt,BlkCt}).PixelIdxList)) = 0;
                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = single(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt});

                % Remove Flagged blobs from regionprops
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(rp_FFT_Filtered_Clean_FlagBlobIdx2_Cell{SkelCt,SegCt,BlkCt}) = [];
                
                % if removed more than 50% of orent bands above, remove all bands
                if (height(find(rp_FFT_Filtered_Clean_FlagBlobIdx2_Cell{SkelCt,SegCt,BlkCt}))/height(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})>0.50)
                    Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = zeros(size(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}),'single');
                    rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = [];
                    MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt} = NaN;
                end

                clearvars clearvars rp_y_range rp_y_range_Temp bCt;

            end
            
            % = Proceed w/ remaining orient bands after clean up
            if ~isempty(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})

                if isnan(MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt}) 
                    MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt} = ...
                        max(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MaxFeretDiameter)./cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MinFeretDiameter),[],'all');
                else
                    MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt} = ...
                        max(MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt},max(cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MaxFeretDiameter)./cat(1,rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}.MinFeretDiameter),[],'all'),'omitmissing');
                end

                if isempty(MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt})
                    MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt} = NaN;
                end

            else

                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = zeros(size(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}),'single');
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = [];
                MaxFeretMaxMinRatio_Cell{SkelCt,SegCt,BlkCt} = NaN;

            end

        end
    end
end
toc


%% Assign weights to slope lines

Weight_Feret_Power = 1;
MaxFeretMaxMinRatio = max(cell2mat(MaxFeretMaxMinRatio_Cell),[],'all','omitmissing');

for SegCt = 1:NumSeg_In % for faster than parfor
    for SkelCt = 1:NumSkel_In
        for BlkCt = 1:NumBlk_In

            if ~isempty(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})

                Img_Temp = zeros(size(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}),'single');
                for bCt = 1:height(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})
                    bWeight_Feret = ...
                        (rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(bCt).MaxFeretDiameter/rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(bCt).MinFeretDiameter)^Weight_Feret_Power;
                    bWeight_Feret = bWeight_Feret/(MaxFeretMaxMinRatio^Weight_Feret_Power);
                    bWeight_FC = mean(Img_FFT_Filtered_Cell_In{SkelCt,SegCt,BlkCt}(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(bCt).PixelIdxList),'all'); 
                    bWeight_pp = abs(mean(Img_Postprocess_Cell_In{SkelCt,SegCt,BlkCt}(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(bCt).PixelIdxList),'all')-1); 
                    bWeight_Orient = width(Img_FFT_Filtered_Cell_In{SkelCt,SegCt,BlkCt})/cosd(abs(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(bCt).Orientation));
                    Img_Temp(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}(bCt).PixelIdxList) = bWeight_Feret*bWeight_FC*bWeight_pp/bWeight_Orient;
                end

                Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = Img_Temp;

            end

        end
    end
end

clearvars Img_Temp bWeight* Img_FFT_Filtered_Cell_Lin Img_FFT_Filtered_Cell_IntensityOutlierHiThresh;


%% Adjust brightness

MaxIntensity_FullCell = max(cellfun(@(x) max(x,[],'all'),Img_FFT_Filtered_Clean_Cell_Out),[],'all');
Img_FFT_Filtered_Clean_Cell_Out = cellfun(@(x) x./MaxIntensity_FullCell,Img_FFT_Filtered_Clean_Cell_Out,'UniformOutput',false);

% = Update regionprops
for SegCt = 1:NumSeg_In 
    for SkelCt = 1:NumSkel_In
        for BlkCt = 1:NumBlk_In

            if ~isempty(rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})
                if ~isreal(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt})
                    1+1
                end

                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = ...
                    regionprops(logical(Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt}),Img_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt},...
                    'Area','PixelIdxList','Orientation','MeanIntensity');
            else
                rp_FFT_Filtered_Clean_Cell_Out{SkelCt,SegCt,BlkCt} = [];
            end

        end
    end
end

