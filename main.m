% This is where we do all the jobs
clear;
clc;

% initialize parameters
alph = 1;
para_d = 1;
para_k = 1;
c_0 = 1;

% set number of elements in x and y, initialize A and b
enumx = 100;
enumy = 20;
a = zeros((enumx + 1) * (enumy + 1));
b = zeros((enumx + 1) * (enumy + 1), 1);

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

% compute base function values and derivatives evaluated
% at each gaussian integration point
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

% iterate over all elements
for i = 1:enumx*enumy
    [nidx, ncoords] = fgetn(i, enumx, enumy);
    a_el = zeros(4, 4);
    [pcoords, det_j, dphi_dx, dphi_dy] = fpdtrans(ncoords, base, db_dxi, db_deta);
    gweit = gweights .* det_j;
    velocity = alph .* (0.25 * 2 ^ 2 - abs(pcoords(:, 2) - 1) .^ 2);

    % the first part, v(y) * (dc / dx)
    for j = 1:4
        for k = 1:4
            a_el(j, k) = a_el(j, k) + base(j, :) .* dphi_dx(k, :) .* velocity' * gweit;
        end
    end

    % second part, D * grad ^ 2 * c
    for j = 1:4
        for k = 1:4
            a_el(j, k) = a_el(j, k) + (para_d .* (dphi_dx(j, :) .* dphi_dx(k, :) + dphi_dy(j, :) .* dphi_dy(k, :))) * gweit;
        end
    end

    % third part k * c
    for j = 1:4
        for k = 1:4
            a_el(j, k) = a_el(j, k) + (para_k .* base(j, :) .* base(k, :)) * gweit;
        end
    end

    % assemble a
    for j = 1:4
        for k = 1:4
            a(nidx(j), nidx(k)) = a(nidx(j), nidx(k)) + a_el(j, k);
        end
    end
end

% set up boundary condition c = c_0 at left
bidxs = fgetbn(enumx, enumy);
for j = 1:length(bidxs)
    a(bidxs(j), :) = 0;
    a(bidxs(j), bidxs(j)) = 1;
    b(bidxs(j)) = c_0;
end

% solve
c = a \ b;

% report result
result = fgenres(c, enumx, enumy);
contourf(result);