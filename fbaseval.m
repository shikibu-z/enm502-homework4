function [base, db_dxi, db_deta] = fbaseval(xi, eta)
% fbaseval - Evaluate base function and derivatives at given
% xi and eta, in UNITE ELEMENT
    % initialize
    base = zeros(4, 1);
    db_dxi = zeros(4, 1);
    db_deta = zeros(4, 1);

    % evaluate base function
    base(1, :) = 0.25 * (1 - xi) * (1 - eta);
    base(2, :) = 0.25 * (1 + xi) * (1 - eta);
    base(3, :) = 0.25 * (1 - xi) * (1 + eta);
    base(4, :) = 0.25 * (1 + xi) * (1 + eta);

    % evaluate derivative of xi
    db_dxi(1, :) = -0.25 * (1 - eta);
    db_dxi(2, :) = 0.25 * (1 - eta);
    db_dxi(3, :) = -0.25 * (1 + eta);
    db_dxi(4, :) = 0.25 * (1 + eta);

    % evaluate derivative of eta
    db_deta(1, :) = -0.25 * (1 - xi);
    db_deta(2, :) = -0.25 * (1 + xi);
    db_deta(3, :) = 0.25 * (1 - xi);
    db_deta(4, :) = 0.25 * (1 + xi);
end