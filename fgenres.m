function result = fgenres(c, enumx, enumy)
% fgenres - Generate value on the node for result report
    result = zeros(enumy + 1, enumx + 1);
    for i = 1:enumy+1
        for j = 1:enumx+1
            nidx = (i - 1) * (enumx + 1) + j;
            result(i, j) = c(nidx);
        end
    end
end