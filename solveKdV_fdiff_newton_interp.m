function [xout, uout] = solveKdV_fdiff_newton_interp(xin, uin, config, gridpoints, c)

% this uses the Newton solver to let us change:
%  N - number of gridpoints
%  c - speed of wave
% large changes in c will likely not work

%% setup

% set up new grid
N = gridpoints;         % number of mesh points, must be even! 
L = xin(end) - xin(1)   % length of domain
h = L/(N-1); 
xout = (0:N-1)'*h;

% interpolate function onto new grid
u = interp1(xin,uin(1:end-1),xout);

% if we don't specify a new c, then take it from the 
% last element of uin
if ~exist('c','var')
    c = uin(end);
end

par.c = c;


%% compute finite difference matrices

[D, D2, D3, D4] = D_fdiff(N, h, config.BC);

%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',50);

% call solve
[uout,fval] = fsolve(@(u) integratedKdV_fdiff(u,D,D2,D3,D4,N,par),u,options);

%% plot results (optional)

% figure;
% plot(xout,uout); % plot solution
% title('Pulse on the full line');

% reappend c to the output vector
uout = [uout ; c];

