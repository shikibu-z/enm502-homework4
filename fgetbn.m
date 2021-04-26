function bidxs = fgetbn(enumx, enumy)
% fgetbn - Get the node index in left boundary.
    bidxs = zeros(enumy + 1, 1);
    for i = 1:enumy+1
        bidxs(i) = 1 + (i - 1) * (enumx + 1);
    end
end