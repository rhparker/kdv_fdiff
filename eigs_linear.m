% eigenvalues of linearization about stationary wave
% uses eigs instead of eig for speed
function [lambda, V, J] = eigs_linear(x, u, config, num, center)
    par.c = u(end);
    udata = u(1:end-1);
    N = length(udata);
    L = x(end) - x(1);             
    h = L/(N-1); 
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
    [F, J] = KdV_fdiff(udata,D,D2,D3,D4,D5,N,par);
    opts.tol = 10^(-10); 
    [V, lambda_D] = eigs(J, num, center, opts);
    lambda = diag(lambda_D);
end
