% This is where we do all the jobs
% set number of elements in x and y first

% initialize gaussian point and weight
gauss_p = [0.7746; 0; -0.7746];
gauss_w = [0.5556; 0.8889; 0.5556];

% initialize bilinear basis functions,
% evaluated at each gaussian point, in UNITE ELEMENT
iter = 0;
base = zeros(4, 9);
db_dxi = zeros(4, 9);
db_deta = zeros(4, 9);

for i = 1:3
    for j = 1:3
        iter = iter + 1;
        [temp_b, temp_dxi, temp_deta] = fbaseval(gauss_p(i), gauss_p(j));
        base(:, iter) = temp_b;
        db_dxi(:, iter) = temp_dxi;
        db_deta(:, iter) = temp_deta;
    end
end