% generate initial data
% takes too long, so we will load it from a file
% N = 1000;                 % number of finite diff gridpts
% iteration = 150;          % number of iterations of continuation code
% [x, uc] = solveKdV_fdiff(N, iterations);

% instead we load initial data from file
% this contains 150 iterations with step 0.25 in c
%  x:  grid of x values, 0 to 50, 1000 grid points
%  uc: continuation data; each column of uc is of the form
%       [u; c], where u is the solution and c is the speed

%% load data and make configuration
% load uc_data;                 % FD, L = 50
load uc_data_100;             % FD, L = 100
% load uc_full;                 % full pulse

% method to use
config.method   = 'fdiff';
config.BC       = 'Neumann';

% specify which index you want in uc
index       = 35;               % L = 50, 100
% index_full  = 48;

% input data to use
u       = uc(:,index);
% u_full  = uc_full(:,index_full);

% ramp up to more grid points
% for increased accuracy
% [xin, uin, c] = solveKdV_fdiff_newton_interp(x, u, config, 2000);
[xin, uin, c] = solveKdV_fdiff_newton_interp(x, u, config, 4000);

% % if we don't do this, then use original values for input
% xin = x;
% uin = u;

% plot(xin,uin(1:end-1));

%% make plot of oscillations
% plot oscillations and get scaled version of u
c = uin(end);
[y, uscaled, start, decay, freq] = osc_plot_c(xin, uin);

%% find the minima and maxima of this
udata = uin(1:end-1);
uscaled = udata.*exp(decay*xin);

% compute finite difference matrices
N = length(xin);
L = xin(end) - xin(1);             
h = L/(N-1); 
D = D_fdiff(N, h, config.BC);
y = D*uscaled./uscaled - decay;

% threshold for a zero of this function
% epsilon = 0.025;               
% epsilon = 0.005;                  % L = 50
epsilon = 0.01;                   % L = 100
zeros = find(abs(y) < epsilon);

% plot(xin, y, '.');
% axis([0 100 -0.5 0.5]);

%% use output of half-wave to make full wave
[x2, u2] = full_wave(xin, uin);

% Newton solver using Neumann BCs
un_out = solveKdV_fdiff_newton(x2, u2, config);

% Newton solver using periodic BCs
p_config.method = 'fdiff'; 
p_config.BC     = 'periodic';
% need to remove last point for periodic BCs
% since it is identified with first point
xper = x2(1:end-1);
uper = [u2(1:end-2); c];
up_out = solveKdV_fdiff_newton(xper, uper, p_config);

% plot(x2, un_out(1:end-1), xper, up_out(1:end-1))

% %% try the double soliton with our zeros
% 
% % for now, do the first 3 zeros
% 
% % for i = 1:4                 % L = 50
% for i = 2:5                 % L = 100
%     x_len   = length(x2) + 1;
%     z       = zeros(i) - 1;
%     ud      = u2(x_len/2 - z + 1:x_len - z);
%     ud_left = flipud(ud);
%     ud_full = [ ud_left(1:end-1) ; ud; c ];
% 
%     ud_out  = solveKdV_fdiff_newton(x2, ud_full, config);
% 
%     figure;
%     plot(x2, ud_full(1:end-1), x2, ud_out(1:end-1));
%     legend('initial guess','Newton solver output')
%     title(strcat('double soliton, speed c =  ',num2str(c)))
% end

%% fun with eigenvalues

x_len   = length(x2) + 1;
z       = zeros(2);
ud      = u2(x_len/2 - z + 1:x_len - z);
ud_left = flipud(ud);
ud_full = [ ud_left(1:end-1) ; ud; c ];

ud_out  = solveKdV_fdiff_newton(x2, ud_full, config);
figure;
plot(x2, ud_full(1:end-1), x2, ud_out(1:end-1));
legend('initial guess','Newton solver output')
title(strcat('double soliton, speed c =  ',num2str(c)))

% don't run this for now, since it takes forever
% and we have the output of eig saved

% [V,J] = eig_linear(x2, ud_full, config)