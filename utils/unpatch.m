function data = unpatch(X, k)
    kernel=k+k+1; 
    data = zeros([kernel, size(X, 1)]);
    tmp = zeros(kernel);
    for i = 1:size(X, 1)
        tmp(:) = X(i, :);
        data(:,:,:,i) = tmp;
    end 
end


