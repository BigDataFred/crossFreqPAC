function [CTFcoord,Gridcoord,VCidx] = convert_MNI2Source_gridCoord(tm,MNIcoord,grid)

CTFcoord = zeros(size(MNIcoord,1),3);
for it = 1:size(MNIcoord,1)
    CTFcoord(it,:) = ft_warp_apply(tm,MNIcoord(it,:),'homogenous');
end;

X = grid.pos(find(grid.inside==1),:);

midx = zeros(size(CTFcoord,1),1);
for it = 1:size(CTFcoord,1)
    [~,midx(it)] = min(sum((X - repmat(CTFcoord(it,:),[size(X,1) 1])).^2,2));
end;

[Gridcoord] = X(midx,:);

[VCidx] =  midx;