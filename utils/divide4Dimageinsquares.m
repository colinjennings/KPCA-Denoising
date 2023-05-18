function [Xcell,icenter,jcenter,kcenter] = divide4Dimageinsquares(S,kernel,mask)

kernel = kernel + (mod(kernel, 2)-1);   % needs to be odd.
N = prod(kernel);
M = size(S,4);
k = (kernel-1)/2; ki = k(1); kj = k(2); kk = k(3); %laterals of the kernel

stats = regionprops(mask, 'BoundingBox');
n = ceil(stats.BoundingBox(4:6) ./ kernel);

icenter = linspace(ceil(stats.BoundingBox(2))+k(2), floor(stats.BoundingBox(2))-k(2) + stats.BoundingBox(5), n(2)); icenter = round(icenter);
jcenter = linspace(ceil(stats.BoundingBox(1))+k(1), floor(stats.BoundingBox(1))-k(1) + stats.BoundingBox(4), n(1)); jcenter = round(jcenter);
kcenter = linspace(ceil(stats.BoundingBox(3))+k(3), floor(stats.BoundingBox(3))-k(3) + stats.BoundingBox(6), n(3)); kcenter = round(kcenter);

[icenter, jcenter, kcenter] = meshgrid(icenter, jcenter, kcenter); icenter = icenter(:); jcenter = jcenter(:); kcenter = kcenter(:);
TotBlock = numel(icenter);
Xcell = cell(TotBlock,1);

for nn = 1:TotBlock
    X = S(icenter(nn)-ki:icenter(nn)+ki, jcenter(nn)-kj:jcenter(nn)+kj, kcenter(nn)-kk:kcenter(nn)+kk, :);
    X = reshape(X, N, M); X = X';
    Xcell{nn} = X;
end
  
end

