% Copyright 2007, David JB Lloyd, Daniele Avitabile, Bjorn Sandstede, Alan R Champneys
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [F,J] = integratedKdV_fdiff(u,D,D2,D3,D4,N,par)
% returns the right-hand side of our equation

%% operator

% 5th order KDV, integrated once
LN = D4 - D2 + sparse(1:N,[1:N],par.c,N,N);
F = LN*u - u.*u;

%% Jacobian
if nargout > 1 
    % 5th order KDV, intergrated once
    J = LN - 2 * sparse(1:N,[1:N],u,N,N);
end

end
