% eigenvalues of linearization about stationary wave
function V = eig_linear(x, u, config)
    par.c = u(end);
    udata = u(1:end-1);
    N = length(udata);
    L = x(end) - x(1);             
    h = L/(N-1); 
    [D, D2, D3, D4] = D_fdiff(N, h, config.BC);
    [F, J] = integratedKdV_fdiff(udata,D,D2,D3,D4,N,par);
    V = eigs(J);
end