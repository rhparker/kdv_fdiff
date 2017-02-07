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
% load uc_data;                   % FD, L = 50
load uc_data_100;             % FD, L = 100
% load uc_full;                 % full pulse

% method to use
config.method   = 'fdiff';
config.BC       = 'Neumann';

% specify which index you want in uc
index = 35;                 % L = 50, c = 1.5650
% index = 35;                 % L = 100, c = 1.5647
% index = 100;                % L = 50, c = 4.1363
% index = 150;                % L = 50, c = 6.3081
% index = 52                  % L = 100, c = 2.1994

% input data to use
u       = uc(:,index);
% u_full  = uc_full(:,index_full);

% ramp up to more grid points for increased accuracy
% standardize c, if needed

% L = 50

% [xin, uin] = solveKdV_fdiff_newton_interp(x, u, config, 2000);
% [xin, uin] = solveKdV_fdiff_newton_interp(x, u, config, 4000);

% L = 100
% here we want to also make c = 1.5650
c = 1.5650
[xin, uin] = solveKdV_fdiff_newton_interp(x, u, config, 2000, c);

% % if we don't do this, then use original values for input
% xin = x;
% uin = u;

plot(xin,uin(1:end-1));

c = uin(end);

%% make plot of oscillations
% plot oscillations and get scaled version of u
c = uin(end);
[y, uscaled, start, decay, freq] = osc_plot_c(xin, uin);

% spacing between pulses is expected to be an integer
% multiple of this
spacing = pi/freq;

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
% epsilon = 0.020;                   % L = 50, c = 6.3081              
% epsilon = 0.008;                   % L = 100, 10000 grid points
% epsilon = 0.003;                   % L = 100, 2000 grid points, c = 1.5650
% epsilon = 0.005;                   % L = 50, c = 1.5650
% epsilon = 0.010;                   % L = 100, c = 1.5650
epsilon = 0.03

% % scatter plot in case we need it to find zeos
% plot(xin, y, '.');
% axis([0 100 -epsilon epsilon]);

minmax = find(abs(y) < epsilon);


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

% which zero to use
% these are spaced by (approximately, since grid) pi/freq
z       = minmax(11);

x_len   = length(x2) + 1;

% right half-wave
ud       = u2(x_len/2 - z + 1:x_len - z);
ud_right = [ud ; c];

% Newton solver on right half wave
ud_right_out  = solveKdV_fdiff_newton(xin, ud_right, config, 1000);
figure;
% plot the half wave
plot(xin, ud_right(1:end-1), xin, ud_right_out(1:end-1));
legend('initial guess','Newton solver output')
title(strcat('double soliton half wave, speed c =  ',num2str(c)))

% plot a zoom-in of half wave
figure;
plot(xin, ud_right(1:end-1), xin, ud_right_out(1:end-1));
legend('initial guess','Newton solver output')
title(strcat('zoom in of double soliton half wave, speed c =  ',num2str(c)))
yzoom = 5e-1;
axis([0 L -yzoom yzoom]);

% construct full wave from half wave
% this is done after the Newton solver is run on the half wave
[~, ud_full] = full_wave(xin, ud_right_out);

% % Newton solver on full wave
% construct full wave from two half waves (before Newton solver done)
% [~,ud_full] = full_wave(xin, ud_right);
% ud_out  = solveKdV_fdiff_newton(x2, ud_full, config, 500);
% figure;
% plot(x2, ud_full(1:end-1), x2, ud_out(1:end-1));
% legend('initial guess','Newton solver output')
% title(strcat('double soliton, full wave, speed c =  ',num2str(c)))

% % Newton solver on full wave constructed from two half-waves
% % construction done after Newton solver run on half wave
% % this is much more accurate
% ud_out  = solveKdV_fdiff_newton(x2, ud_full_fromhalf, config, 500);
% figure;
% plot(x2, ud_full_fromhalf(1:end-1), x2, ud_out(1:end-1));
% legend('initial guess','Newton solver output')
% title(strcat('double soliton, full wave, speed c =  ',num2str(c)))

% don't run this for now, since it takes forever
% and we have the output of eig saved

% [V,J] = eig_linear(x2, ud_out, config)