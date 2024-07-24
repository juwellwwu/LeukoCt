function D2 = fMeasureChessboardDist(ZI,ZJ)

nx = 1;
px = size(ZI,2); 

ny = size(ZJ,1);  
py = px;
outClass = class(ZI);

dsq_2D = zeros(ny,px,outClass);
for q = 1:px 
    dsq_2D(:,q) = abs(repmat(ZI(1,q),ny,1)-ZJ(:,q)); 
end

D2 = max(dsq_2D,[],2);

end