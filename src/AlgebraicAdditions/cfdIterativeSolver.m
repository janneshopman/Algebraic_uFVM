function [x, relres, iter, resvec] = iterSol(A, M, b, tol, maxIter, x)
    nc = length(b);

    if ~exist('tol', 'var')
        tol = 1E-6;
    end
    
    if ~exist('maxIter', 'var')
        maxIter = 100;
    end

    if ~exist('x', 'var')
        x = zeros(nc, 1);
    end

    resvec = zeros(maxIter + 1, 1);

    normb = vecnorm(b);

    r = b - A*x;
    normr = vecnorm(r);

    resvec(1) = normr;
    relres = normr/normb;

    iter = 0;

    while(relres > tol && iter < maxIter)
        iter = iter + 1;

        x = x + M\r;

        r = b - A*x;
        normr = vecnorm(r);

        resvec(iter + 1) = normr;
        relres = normr/normb;
    end

    resvec = nonzeros(resvec);
end