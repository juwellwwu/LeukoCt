function [] = fSaveImgDRXCellMontage(SaveFilePath_Input,Img_Cell1,varargin)

fprintf('Saving Montage of xDyT Image Cells ...\n');

NumSeg_Input = size(Img_Cell1,2);

for WhichSeg=1:NumSeg_Input

    Target_ImgSize_C1 = max(cell2mat(cellfun(@size,Img_Cell1(:,WhichSeg,:),'UniformOutput',false)),[],1);
    Target_ImgSize_C1 = max(Target_ImgSize_C1,[],3);
    padarray_wrapper_C1 = @(x) padarray(x,Target_ImgSize_C1-size(x),0,'post');
    Img_Cell1_Montage = ...
        cellfun(padarray_wrapper_C1,Img_Cell1(:,WhichSeg,:),'UniformOutput',false);
    Img_Cell1_Montage = cat(2,Img_Cell1_Montage{:});
    Img_Cell1_Montage(:,cumsum(repmat(Target_ImgSize_C1(1,2),1,size(Img_Cell1,1)*size(Img_Cell1,3)-1))) = 0;

    if width(varargin)>0 
        Img_Cell2 = varargin{1,1};
        Target_ImgSize_C2 = max(cell2mat(cellfun(@size,Img_Cell2(:,WhichSeg,:),'UniformOutput',false)),[],1);
        Target_ImgSize_C2 = max(Target_ImgSize_C2,[],3);
        if sum(Target_ImgSize_C2)>0
            padarray_wrapper_C2 = @(x) padarray(x,Target_ImgSize_C2-size(x),0,'post');
            Img_Cell2_Montage = ...
                cellfun(padarray_wrapper_C2,Img_Cell2(:,WhichSeg,:),'UniformOutput',false);
            Img_Cell2_Montage = cat(2,Img_Cell2_Montage{:});
            Img_Cell2_Montage(:,cumsum(repmat(Target_ImgSize_C2(1,2),1,size(Img_Cell2,1)-1))) = 0;
        else
            Img_Cell2_Montage = zeros(1,size(Img_Cell1_Montage,2),'single');
        end
    else
        Img_Cell2_Montage = zeros(1,size(Img_Cell1_Montage,2),'single');
    end

    if width(varargin)>1 
        Img_Cell3 = varargin{1,2};
        Target_ImgSize_C3 = max(cell2mat(cellfun(@size,Img_Cell3(:,WhichSeg,:),'UniformOutput',false)),[],1);
        Target_ImgSize_C3 = max(Target_ImgSize_C3,[],3);
        if sum(Target_ImgSize_C3)>0
            padarray_wrapper_C3 = @(x) padarray(x,Target_ImgSize_C3-size(x),0,'post');
            Img_Cell3_Montage = ...
                cellfun(padarray_wrapper_C3,Img_Cell3(:,WhichSeg,:),'UniformOutput',false);
            Img_Cell3_Montage = cat(2,Img_Cell3_Montage{:});
            Img_Cell3_Montage(:,cumsum(repmat(Target_ImgSize_C3(1,2),1,size(Img_Cell3,1)-1))) = 0;
        else
            Img_Cell3_Montage = zeros(1,size(Img_Cell1_Montage,2),'single');
        end
    else
        Img_Cell3_Montage = zeros(1,size(Img_Cell1_Montage,2),'single');
    end

    Img_Montage = zeros((size(Img_Cell1_Montage,1)+size(Img_Cell2_Montage,1)+size(Img_Cell3_Montage,1)),...
        max([size(Img_Cell1_Montage,2);size(Img_Cell2_Montage,2);size(Img_Cell3_Montage,2)],[],1),...
        'single');
    Img_Montage(1:size(Img_Cell1_Montage,1),1:size(Img_Cell1_Montage,2)) = Img_Cell1_Montage;
    Img_Montage((size(Img_Cell1_Montage,1)+1):(size(Img_Cell1_Montage,1)+size(Img_Cell2_Montage,1)),1:size(Img_Cell2_Montage,2)) = ...
        Img_Cell2_Montage;
    Img_Montage((size(Img_Cell1_Montage,1)+size(Img_Cell2_Montage,1)+1):end,1:size(Img_Cell3_Montage,2)) = ...
        Img_Cell3_Montage;
    Img_Montage(size(Img_Cell1_Montage,1)+size(Img_Cell2_Montage,1)-1:size(Img_Cell1_Montage,1)+size(Img_Cell2_Montage,1)+1,:) = 0;

    %         figure;
    %         imshow(Img_Montage,'InitialMagnification',100);

    imwrite(Img_Montage,...
        strcat(SaveFilePath_Input,'Img_xDyT_Montage_TimeSeg',num2str(WhichSeg,'%06d'),'.jpg'));

end