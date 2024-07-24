function [FrIntensitySum_Struct_Out] = fCreateFrIntensitySumStruct(ImgStack_In,Num_Blk_or_WinLength_In)

% = Define structure 
FrIntensitySum_Struct_Out = struct;

FrIntensitySum_Struct_Out.IPZmax =  -9999*ones(size(ImgStack_In,3),Num_Blk_or_WinLength_In); 
FrIntensitySum_Struct_Out.IPRaw =  -9999*ones(size(ImgStack_In,3),Num_Blk_or_WinLength_In); 
FrIntensitySum_Struct_Out.IPBkgd =  -9999*ones(size(ImgStack_In,3),Num_Blk_or_WinLength_In);
FrIntensitySum_Struct_Out.IPBkgdRmv =  -9999*ones(size(ImgStack_In,3),Num_Blk_or_WinLength_In);
FrIntensitySum_Struct_Out.IPMAD =  -9999*ones(size(ImgStack_In,3),Num_Blk_or_WinLength_In);

FrIntensitySum_Struct_Out.Img_Mask_Cell = cell(1,Num_Blk_or_WinLength_In); 

FrIntensitySum_Struct_Out.Img_Mask_Overlay_Cell =... 
    repmat({zeros(height(ImgStack_In),width(ImgStack_In),3,'uint8')},1,Num_Blk_or_WinLength_In); 

FrIntensitySum_Struct_Out.MaskSkelLength_px = zeros(1,Num_Blk_or_WinLength_In); 
FrIntensitySum_Struct_Out.MaskSkelLength_um = zeros(1,Num_Blk_or_WinLength_In);

FrIntensitySum_Struct_Out.MaskWBCAreaFrac = zeros(1,Num_Blk_or_WinLength_In);

FrIntensitySum_Struct_Out.VesselDiameter_px = zeros(1,Num_Blk_or_WinLength_In);
FrIntensitySum_Struct_Out.VesselDiameter_um = zeros(1,Num_Blk_or_WinLength_In);

FrIntensitySum_Struct_Out.FocusDiameterIdxMean = zeros(1,Num_Blk_or_WinLength_In);

FrIntensitySum_Struct_Out.CellCtRaw = struct;
FrIntensitySum_Struct_Out.CellCtRaw.CellCt = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCtRaw.Method = 'None';
FrIntensitySum_Struct_Out.CellCtRaw.Optimized = -9999*ones(1,Num_Blk_or_WinLength_In,'single'); 
FrIntensitySum_Struct_Out.CellCtRaw.ZScore = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCtRaw.pkParam = repmat({zeros(0,0,'single')},1,Num_Blk_or_WinLength_In); 
FrIntensitySum_Struct_Out.CellCtRaw.rpGateXY = repmat({zeros(0,0,'single')},1,Num_Blk_or_WinLength_In); 
FrIntensitySum_Struct_Out.CellCtRaw.IPCC = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCtRaw.IPMADCC = -9999*ones(1,Num_Blk_or_WinLength_In,'single'); 
FrIntensitySum_Struct_Out.CellCtRaw.IPBkgdRmvPkAg = -9999*ones(1,Num_Blk_or_WinLength_In,'single'); 
FrIntensitySum_Struct_Out.CellCtRaw.BrightRankVal = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCtRaw.BrightRank = -9999*ones(1,Num_Blk_or_WinLength_In,'single');

FrIntensitySum_Struct_Out.CellCt = struct;
FrIntensitySum_Struct_Out.CellCt.CellCt = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCt.Method = 'None';
FrIntensitySum_Struct_Out.CellCt.Optimized = -9999*ones(1,Num_Blk_or_WinLength_In,'single'); 
FrIntensitySum_Struct_Out.CellCt.ZScore = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCt.pkParam = repmat({zeros(0,0,'single')},1,Num_Blk_or_WinLength_In); 
FrIntensitySum_Struct_Out.CellCt.rpGateXY = repmat({zeros(0,0,'single')},1,Num_Blk_or_WinLength_In); 
FrIntensitySum_Struct_Out.CellCt.IPCC = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCt.IPMADCC = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCt.IPBkgdRmvPkAg = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCtRaw.BrightRankVal = -9999*ones(1,Num_Blk_or_WinLength_In,'single');
FrIntensitySum_Struct_Out.CellCtRaw.BrightRank = -9999*ones(1,Num_Blk_or_WinLength_In,'single');