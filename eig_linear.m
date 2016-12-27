% eigenvalues of linearization about stationary wave
function [V, J] = eig_linear(x, u, config)
    par.c = u(end);
    udata = u(1:end-1);
    N = length(udata);
    L = x(end) - x(1);             
    h = L/(N-1); 
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
    [F, J] = KdV_fdiff(udata,D,D2,D3,D4,D5,N,par);
    opts.tol = 10^(-10); 
    V = eig(full(J));
%     V = eigs(J, 10, 0, opts);
end