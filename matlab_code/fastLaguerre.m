function L = fastLaguerre(n, alpha, x)
    % (just to do the work faster than built in laguerre matlab
    % Calculate the generalized Laguerre polynomial of degree n and parameter alpha
    L = zeros(size(x));
    for k = 0:n
        coeff = (-1)^k * factorial(n + alpha) / (factorial(k) * factorial(n - k) * factorial(alpha + k));
        L = L + coeff * x.^k;
    end
end