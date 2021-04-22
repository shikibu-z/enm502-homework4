function [nidxs, ncoords] = fgetn(eidx, enumx, enumy)
% fgetn - Get the nodes indices and coordinates of an element.
% This function will return the four corner nodes indices and
% coordinates of a given element index. Assuming the domain is
% a rectangle of 10 by 2 and left corner placed at the origin.
% The elements are divided evenly on the domain as rectangles.
    % initialize
    nidxs = [];
    ncoords = [];

    % compute element position
    ecol = mod(eidx, enumx);
    if ecol == 0
        ecol = enumx;
        erow = fix(eidx / enumx);
    else
        erow = fix(eidx / enumx) + 1;
    end

    % calculate node indices
    nidx_base = (erow - 1) * (enumx + 1) + ecol;
    nidxs(1, :) = nidx_base;
    nidxs(2, :) = nidx_base + 1;
    nidxs(3, :) = nidx_base + (enumx + 1);
    nidxs(4, :) = nidx_base + (enumx + 1) + 1;

    % calculate node coordinates
    nbase_coor = [(ecol - 1) * 10 / enumx, (erow - 1) * 2 / enumy];
    ncoords(1, :) = nbase_coor;
    ncoords(2, :) = [nbase_coor(1) + (10 / enumx), nbase_coor(2)];
    ncoords(3, :) = [nbase_coor(1), nbase_coor(2) + (2 / enumy)];
    ncoords(4, :) = [nbase_coor(1) + (10 / enumx), nbase_coor(2) + (2 / enumy)];
end