function [uout, c] = solveKdV_fdiff_newton(xold,uold)

% - Solves 1D quadratic-cubic Swift-Hohenberg equation
% - Finds localised pulse in snaking region and computes its stability
% - Uses finite differences and sparse matrices
% - Requires optimization toolbox and external routine SH_rhs_finite_differences.m
%
% Copyright 2007, David JB Lloyd, Daniele Avitabile, Bjorn Sandstede, Alan R Champneys
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% D N
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% setup

% here use the Newton solver on the same domain as given
N = length(xold);
L = xold(end) - xold(1);             
h = L/(N-1); 

par.c = uold(end);
c = par.c;

u = uold(1:end-1);

%% compute finite difference matrices

% d_x
D = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N); 
D = (D - D');
% Neumann boundary conditions (first derivative 0 at L and R)
D(1,2) = 0; 
D(N,N-1) = 0;
D = D./(2 * h);

% d_xx
D2 = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) - sparse(1:N,[1:N],ones(N,1),N,N);
D2 = (D2 + D2');
% Neumann boundary conditions (inherited from D)
D2(1,2) = 2; 
D2(N,N-1) = 2;
D2 = D2/h^2;

% d_xxx
D3 = -2 * sparse(1:N-1,[2:N],ones(N-1,1),N,N) + sparse(1:N-2,[3:N],ones(N-2,1),N,N);
D3 = (D3 - D3');
% Neumann boundary conditions (third derivative 0 at L and R)
D3(1,2)   = 0; D3(1,3)   = 0;
D3(N,N-1) = 0; D3(N,N-2) = 0;
% Neumann boundary conditions (inherited from D)
D3(2,2)= -1; D3(N-1,N-1) = 1;
D3 = D3/(2 * h^3);

% d_xxxx
D4 = sparse(1:N,[1:N],3 * ones(N,1),N,N) - sparse(1:N-1,[2:N],4 * ones(N-1,1),N,N) + sparse(1:N-2,[3:N],ones(N-2,1),N,N);
D4 = (D4 + D4');
% Neumann boundary conditions (inherited from above)
D4(1,2) = -8; D4(1,3) = 2;
D4(2,2) = 7;
D4(N,N-1) = -8; D4(N,N-2) = 2;
D4(N-1,N-1) = 7;
D4 = D4/h^4;

% identity matrix
e = sparse(1:N,[1:N],1,N,N);

%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',50);

% call solve
[uout,fval] = fsolve(@(u) integratedKdV_fdiff(u,D,D2,D3,D4,N,par),u,options);

%% plot results (optional)

% figure;
% plot(x,uout); % plot solution
% title('Pulse on the full line');

