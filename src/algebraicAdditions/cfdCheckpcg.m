function [dtp, relres, iter, resvec] = cfdCheckpcg(A, b, tol, maxIter, preL, preR, x)

[dtp, flag, relres, iter, resvec] = pcg(A, b, tol, maxIter, preL, preR, x);

if flag == 0
    fprintf('\npcg: relres: %3.2e, iter: %d\n', relres, iter);
end

switch flag
    case 1
        error('PCG iterated MAXIT times but did not converge.');
    case 2
        error('preconditioner M was ill-conditioned.');
    case 3 
        warning('PCG stagnated (two consecutive iterates were the same).');
    case 4 
        warning('one of the scalar quantities calculated during PCG became too small or too large to continue computing.');
    otherwise
end