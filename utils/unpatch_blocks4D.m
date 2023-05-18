function [ Signal ] = unpatch_blocks4D(Xcell,icenter,jcenter,kcenter, kernel,size_)

TotBlocks = numel(icenter(:));
kernel = kernel + (mod(kernel, 2)-1);   % needs to be odd.
k = (kernel-1)/2; ki = k(1); kj = k(2); kk = k(3); %laterals of the kernel
Signal = zeros(size_);

for nn = 1:TotBlocks
     Signal(icenter(nn)-ki:icenter(nn)+ki,jcenter(nn)-kj:jcenter(nn)+kj,kcenter(nn)-kk:kcenter(nn)+kk, :) = unpatch(Xcell{nn}, k);
end

end
