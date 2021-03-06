% eigenvalues of linearization about stationary wave
function [lambda, V, J] = eig_linear(x, u, config)
    par.c = u(end);
    udata = u(1:end-1);
    N = length(udata);
    L = x(end) - x(1);             
    h = L/(N-1); 
    [D, D2, D3, D4, D5] = D_fdiff(N, h, config.BC);
    [F, J] = KdV_fdiff(udata,D,D2,D3,D4,D5,N,par);
    opts.tol = 10^(-10); 
    [V, LD] = eig(full(J));
%     [V, LD] = eig(full(J),'nobalance');
    lambda = diag(LD);
end
