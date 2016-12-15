function [y,uscaled, start] = osc_plot_c(x, u, c)
% find the spatial eigenvalues
% roots of nu^4 - nu^2 + c == 0
nu = roots([1 0 -1 0 c]);

% oscillations frequency is imag(nu)
% oscillations decay with constant re(nu)
decay = abs(real(nu(1)));
freq  = abs(imag(nu(1)));

% where to start plot
start = floor(length(x)/200);
y=x(start:end);

% scale solution by exp(decay) to recover oscillations
uscaled=u(start:end).*exp(decay*y);
umax=max(uscaled);

% plot along with sine function of same scaling
figure;
plot(y,uscaled,y,umax*sin(y*freq)); 

legend('rescaled solution','sine function')
title(strcat('scale to see oscillations, speed c =  ',num2str(c)))
end
