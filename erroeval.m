% initialize results
l2_norms = zeros(21, 1);
esize = zeros(21, 1);
iter = 0;

% iterate over different gird size
for i = 20:5:120
    iter = iter + 1;

    % original
    enumx_old = i;
    enumy_old = i / 5;
    result_old = main(13, 1, 1, 3, enumx_old, enumy_old);

    % increased mesh size
    enumx_new = enumx_old * 2;
    enumy_new = enumy_old * 2;
    result_new = main(13, 1, 1, 3, enumx_new, enumy_new);

    % get the value on node point from original position
    result_acc = zeros(enumy_old + 1, enumx_old + 1);
    for j = 1:(enumy_old+1)
        for k = 1:(enumx_old+1)
            result_acc(j, k) = result_new((j - 1) * 2 + 1, (k - 1) * 2 + 1);
        end
    end

    % report values
    difference = sum((result_acc - result_old) .^ 2, 'all');
    l2_norms(iter) = sqrt(difference) / ((enumx_old + 1) * (enumy_old + 1));
    esize(iter) = (10 / enumx_old) * (2 / enumy_old);
end