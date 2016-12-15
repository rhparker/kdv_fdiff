function osc(x, uc, index)
%% make plot of oscillations
% extract initial data
uin = uc(1:end,index);

% ramp up to 5000 grid points
[xfine, ufine, c] = solveKdV_fdiff_newton_interp(x,uin);
[y,uscaled, start] = osc_plot_c(xfine, ufine, c);

u2 = vertcat( zeros(start-1,1), uscaled )

%% find locations of minima and maxima of oscillations
[peaks, max_locs] = findpeaks(u2);
[peaks, min_locs] = findpeaks(-u2);
min1 = min_locs(1);
max1 = max_locs(1);

% piece together two halves to get whole wave
uleft = flipud(ufine);
ufull = [ uleft(1:end-1) ; ufine ];

% piece together domain as well, for plotting purposes
xleft = -flipud(xfine);
xfull = [ xleft(1:end-1); xfine ];

% test this on the Newton solver
ufull_newton = solveKdV_fdiff_newton(xfull, [ufull;c])
figure;
plot(xfull,ufull,xfull,ufull_newton);

%% construct a double soliton
% join at the first maximum and the first minimum
% of the oscillations; since only plotting half the wave,
% we are technically not joining anything. Append speed c
% for Newton solver
ud_min = [ufull(5001 - min1:10000 - min1); c];
ud_max = [ufull(5001 - max1:10000 - max1); c];

% run both of these through the Newton solver
[xout,ud_min_out] = solveKdV_fdiff_newton_interp(xfine,ud_min);
[xout,ud_max_out] = solveKdV_fdiff_newton_interp(xfine,ud_max);

% plot results
figure;
plot(xfine, ud_min(1:end-1), xfine, ud_min_out);
legend('initial guess (first min)','result from Newton solver')
title('double soliton experiment, join at first min');

figure;
plot(xfine, ud_max(1:end-1), xfine, ud_max_out);
legend('initial guess (first max)','result from Newton solver')
title('double soliton experiment, join at first max');

