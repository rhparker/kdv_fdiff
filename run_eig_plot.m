%% plots eigenvalues of linearization about various waves
% loads eigenvalues from file since
% eig takes forever to run

% ranges for plots
% standard ranges
xrange      = 1e-01;
yrange      = 1e9;
% range to see distinct eigenvalues near real axis
yrange_zoom = 1e-4;

% single soliton, for comparison purposes
% load eig5Ns;

% figure;
% plot(x2, u2(1:end-1));
% title(strcat('Single soliton soliton, speed c =  ',num2str(c)) )
% 
% figure
% plot(V,'.');
% title('Eigenvalues for single soliton solution');
% axis([-xrange xrange -yrange yrange]);

% double soliton
% % load eig5Nd1;
% load eig5Nd1_100;
% load eig5Nd2;
% load eig5Nd3;
% load eig5Nd4;
load eigL100_d10;

% figure;
% plot(x2, ud_full(1:end-1));
% title(strcat('Double soliton soliton, speed c =  ',num2str(c)) )

figure;
plot(lambda,'.');
title('Eigenvalues for double pulse, 1st min/max (L = 100, c = 1.5650)');
% axis([-xrange xrange -yrange yrange]);

% figure;
% plot(lambda,'.');
% title('Eigenvalues for double soliton solution, zoomed in');
% axis([-xrange xrange -yrange_zoom yrange_zoom]);
% 
% % use for real eigenfunctions
% figure;
% plot(x2,real(V))
% legendcell = cellstr(num2str(lambdaV'))
% legend(legendcell);
% title('Eigenvectors for double soliton solution');

% % use for complex eigenfunctions
% figure;
% plot(x2,real(V))
% legendcell = cellstr(num2str(lambdaV))
% legend(legendcell);
% title('Real part of eigenvectors for double soliton solution');
% 
% figure;
% plot(x2,imag(V))
% legend(legendcell);
% title('Imaginary part of eigenvectors for double soliton solution');
% 
% figure;
% plot(x2,abs(V))
% legend(legendcell);
% title('Absolute value of eigenvectors for double soliton solution');

