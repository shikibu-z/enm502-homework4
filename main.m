% This is where we do all the jobs
% set number of elements in x and y first
enumx = 5;
enumy = 3;
a = zeros((enumx + 1) * (enumy + 1));

% initialize gaussian point and weight
gauss_p = [0.7746; 0; -0.7746];
gauss_w = [0.5556; 0.8889; 0.5556];

% initialize bilinear basis functions,
% evaluated at each gaussian point, in UNITE ELEMENT
iter = 0;
base = zeros(4, 9);
db_dxi = zeros(4, 9);
db_deta = zeros(4, 9);
gweights = zeros(9, 1);

for i = 1:3
    for j = 1:3
        iter = iter + 1;
        [temp_b, temp_dxi, temp_deta] = fbaseval(gauss_p(i), gauss_p(j));
        base(:, iter) = temp_b;
        db_dxi(:, iter) = temp_dxi;
        db_deta(:, iter) = temp_deta;
        gweights(iter) = gauss_w(i) * gauss_w(j);
    end
end

for i = 1:enumx*enumy
    [nidx, ncoords] = fgetn(i, enumx, enumy);
    a_el = zeros(4, 4);
    [pcoords, det_j, dphi_dx, dphi_dy] = fpdtrans(ncoords, base, db_dxi, db_deta);
    gweit = gweights .* det_j;
    iter = 0;
    for j = 1:4
        for k = 1:4
            % this is only for grad^2 c
            a_el(j, k) = (dphi_dx(j, :) .* dphi_dx(k, :) + dphi_dy(j, :) .* dphi_dy(k, :)) * gweit;
        end
    end
    for j = 1:4
        for k = 1:4
            a(nidx(j), nidx(k)) = a(nidx(j), nidx(k)) + a_el(j, k);
        end
    end
end