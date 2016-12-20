function osc(x, uc, index, config)
%% make plot of oscillations
% extract initial data
uin = uc(1:end,index);

% ramp up to 5000 grid points
[xfine, ufine, c] = solveKdV_fdiff_newton_interp(x,uin, config);

% plot oscillations and get scaled version of u
[y,uscaled, start] = osc_plot_c(xfine, ufine, c);

% this scaled version has same domain as ufine
uscaled2 = [ zeros(start-1,1); uscaled ];

%% find locations of minima and maxima of oscillations
[peaks, max_locs] = findpeaks(uscaled2);
[peaks, min_locs] = findpeaks(-uscaled2);
min1 = min_locs(1);
max1 = max_locs(1);

% piece together two halves to get whole wave
uleft = flipud(ufine);
ufull = [ uleft(1:end-1) ; ufine ];

% piece together domain as well
xleft = -flipud(xfine);
xfull = [ xleft(1:end-1); xfine ];

% test to make sure full domain works on the Newton solver
ufull_newton = solveKdV_fdiff_newton(xfull, [ufull;c], config);
figure;
plot(xfull,ufull,xfull,ufull_newton);
legend('constructed solution','result from Newton solver')
title('solution on full domain');

%% construct a double soliton
% join at the first maximum and the first minimum
% of the oscillations; since only plotting half the wave,
% we are technically not joining anything. Append speed c
% for Newton solver
ud_min = ufull(5001 - min1:10000 - min1);
ud_max = ufull(5001 - max1:10000 - max1);

ud_min_left = flipud(ud_min);
ud_min_full = [ ud_min_left(1:end-1) ; ud_min ];

ud_max_left = flipud(ud_max);
ud_max_full = [ ud_max_left(1:end-1) ; ud_max ];

%% run both of these through the Newton solver
% % half versions
% ud_min_out = solveKdV_fdiff_newton(xfine,[ud_min;c]);
% ud_max_out = solveKdV_fdiff_newton(xfine,[ud_max;c]);
% xplot = xfine;
% ud_min_start = ud_min;
% ud_max_start = ud_max;

% full versions
ud_min_out = solveKdV_fdiff_newton(xfull,[ud_min_full;c], config);
ud_max_out = solveKdV_fdiff_newton(xfull,[ud_max_full;c], config);
xplot = xfull;
ud_min_start = ud_min_full;
ud_max_start = ud_max_full;

% plot results
figure;
plot(xplot, ud_min_start, xplot, ud_min_out);
legend('initial guess (first min)','result from Newton solver')
title('double soliton experiment, join at first min');

figure;
plot(xplot, ud_max_start, xplot, ud_max_out);
legend('initial guess (first max)','result from Newton solver')
title('double soliton experiment, join at first max');

