function [x,contdata] = solveKdV_fdiff(gridpts, iterations)
%% setup

N = gridpts;
L = 50;             % domain truncation
h = L/(N-1); 
x = (0:N-1)'*h;
% x = (-(N-1):N-1)'*h
% N = 2*N - 1;

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
% Enforce second derivative 0 on both sides
% D2(1,1) = 0; D2(1,2)   = 0;
% D2(N,N) = 0; D2(N,N-1) = 0;
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
% Enforce fourth derivative 0 on both sides
% D4(1,1) = 0; D4(1,2) = 0;   D4(1,3) = 0;
% D4(N,N) = 0; D4(N,N-1) = 0; D4(N,N-2) = 0;
D4 = D4/h^4;

% identity matrix
e = sparse(1:N,[1:N],1,N,N);


%% initial conditions

par.c = 36/169;

% 5th order KDV solution, for c = 36/169
u = (105/338)*sech( x / (2 * sqrt(13) ) ).^4;

%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',50);

% call solve
[uout,fval] = fsolve(@(u) integratedKdV_fdiff(u,D,D2,D3,D4,N,par),u,options);

%% plot results

% figure;
% plot(x,uout); % plot solution
% title('Pulse on the full line');

%% secant continuation of solution to 5th order KdV

%% continuation parameters
% value to start at
par.c = 36/169;

% appears that 0.25 is the biggest step we can takef
contPar.numContSteps = iterations;
contPar.ds = 2.5e-1;         % continuation step size: should always be positive!
contPar.initial_ds = 2.5e-1; % initial step: sign determines direction of continuation
contPar.Name = 'c';        % continuation parameter continued in

%% system parameters

% initial condition is output from above
u0 = uout

%% find two initial points on the bifurcation curve

% first point
options = optimset('Display','iter','Algorithm','levenberg-marquardt','MaxIter',30,'Jacobian','on');
[u1,fval,exitflag,output,jacobian1]  = fsolve(@(u) integratedKdV_fdiff(u,D,D2,D3,D4,N,par),u0,options);

v0 = [u1; getfield(par,contPar.Name)];  % v0 is the first point with parameter name

parameter = getfield(par,contPar.Name); % Set a few facility vectors
normu     = norm(u1);
contdata  = v0;

% second point
par = setfield(par,contPar.Name,getfield(par,contPar.Name)+contPar.initial_ds); % increase par by ds
[u2,fval,exitflag,output,jacobian1]  = fsolve(@(u) integratedKdV_fdiff(u,D,D2,D3,D4,N,par),u1,options);

v1 = [u2; getfield(par,contPar.Name)]; % v1 is the second point with parameter name

parameter = [parameter getfield(par,contPar.Name)];
normu     = [normu     norm(u2)];
contdata  = [contdata  v1];


%% Continuation code
% At each continuation step
for i = 1:contPar.numContSteps

  % Predictor
  v = v1 + (v1 - v0)/norm(v1 - v0, 2) * contPar.ds;
  
  disp(['Predictor = ',num2str(v(end))]);
  % Call fsolve 
  [v,res,exitflag,output,jacobian2] = fsolve(@(v) FixedPointSecantPredictorCorrector_fdiff(v,v1,v0,L,D,D2,D3,D4,N,par,contPar),v,options); 
  
  disp(['Step = ',int2str(i),' Parameter = ',num2str(v(end)), ' Norm(residual,inf) = ',num2str(norm(res,inf))]);
  % Update output parameters
  parameter = [parameter v(end)];
%  normu     = [normu norm(v(1:end-1))];
  normu     = [normu norm(v(1:end-1))];
  contdata  = [contdata v];

  plot(parameter,normu,'-o');   % plot the bifurcation diagram on the fly
  xlabel(contPar.Name);
  ylabel('Norm of the solution');drawnow;
  
  % Prepare for the next continuation step
  v0 = v1;
  v1 = v;
  
end
