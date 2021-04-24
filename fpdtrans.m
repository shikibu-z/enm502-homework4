function [pcoords, det_j, dphi_dx, dphi_dy] = fpdtrans(ncoords, base, db_dxi, db_deta)
% fpdtrans - Calculate partial deratives of transformation
% at 9 gaussian quadrature points
    % initialize
    pcoords = zeros(9, 2);
    det_j = zeros(9, 1);
    dphi_dx = zeros(4, 9);
    dphi_dy = zeros(4, 9);

    % iterate 9 gaussian points
    for i = 1:9
        % transform unite element xi and eta into global coordinates
        pcoords(i, :) = [sum(ncoords(:, 1) .* base(:, i)), sum(ncoords(:, 2) .* base(:, i))];
        % compute pds for coordinate transform
        dx_dxi = sum(ncoords(:, 1) .* db_dxi(:, i));
        dx_deta = sum(ncoords(:, 1) .* db_deta(:, i));
        dy_dxi = sum(ncoords(:, 2) .* db_dxi(:, i));
        dy_deta = sum(ncoords(:, 2) .* db_deta(:, i));
        det_val = dx_dxi * dy_deta - dx_deta * dy_dxi;
        % compute determinate J for each xi and eta
        det_j(i) = abs(det_val);
        % compute dphi_dx and dphi_dy in terms of xi and eta
        dphi_dx(:, i) = (db_dxi(:, i) .* dy_deta - db_deta(:, i) .* dy_dxi) / det_val;
        dphi_dy(:, i) = (db_deta(:, i) .* dx_dxi - db_dxi(:, i) .* dx_deta) / det_val;
    end
end