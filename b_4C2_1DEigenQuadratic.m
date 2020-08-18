%===========================================================
% Example 4C.2
% Atalla
% Solve a 1D Acoustic Eigenvalue Problem using quadratic elements
% Boundary conditions: Fixed pressure at both ends%
%
%===========================================================
clear all; close all; clc
%
% Problem data
%
Lx=1; % Length of the domain; fixed arbitrary to 1
c0=342; % Speed of sound (m/s)
%
% Step 1: Mesh
%
ne= input(' Number of quadratic element: '); % number of elements
nnt = 2*ne+1; % Total number of nodes
h=Lx/ne; % Length of the elements
x = (0:h/2:Lx); % Coordinates table
%
% Step 2: Compute Elementary matrices
%
Ke=c0^2*[7,-8,1;-8,16,-8;1,-8,7]/(3*h);
Me=[4,2,-1;2,16,2;-1,2,4]*h/30;
%
% Step 3: Assembling
%
I=eye(3,3);
K=zeros(nnt,nnt); M = zeros(nnt,nnt);
for ie=1: ne
 L=zeros(3,nnt); L(:,2*ie-1:2*ie+1)=I; % Location matrix for element ie
 K = K + L'*Ke*L;
 M = M + L'*Me*L;
end
%
% Step 4: Boundary conditions: p(1) = p(nnt) = 0
%
K = K(2:nnt-1,2:nnt-1);
M = M(2:nnt-1,2:nnt-1);
ndof = nnt-2; % Final number of unknown (equations)
%
% Step 5: Compute the eigenvalues and eigenvectors
%
[V,D]=eig(K,M); D= sqrt(D);
%
% Normalisation et ordre des vecteurs propres
%
Nor=V'*M*V; V=V*sqrt(inv(Nor)); % Normalization of the eigenvectors
for idof = 1: ndof
    w(idof) = D(idof,idof);
end
[w,ind]=sort(w); % Sort the eigenvalues in ascending order
%
% Step 6: Comparison with the exact solution for a selected number of modes
%
%
nm = min(4, ndof); % Here we select the first modes up to the fourth
theo = c0*[1:ndof]*pi/Lx;
Err = abs(theo-w)./theo;
fprintf('\n ************Results for %d quadratic elements************', ne);
fprintf('\n ================================================');
fprintf('\n Mode FEM (Hz) Exact (Hz) Error(%%)');
fprintf('\n ----------------------------------------------');
for m=1:nm
 fprintf('\n %d %7.2f %7.2f %3.1f ', m, w(m)/2/pi,theo(m)/2/pi, Err(m)*100);
end
fprintf('\n ================================================\n');
fprintf('\n');
%
% comparison of the mode shapes
%
xt=(0:0.025:1) * Lx; % Mesh to visualize the mode shapes
for m=1:nm
     figure;
     yt= sqrt(2)*sin(m*pi*xt); plot(xt,yt), hold on, % Theoretical solution for  
     y=[0,real(V(:,ind(m))'),0]; 
     if (yt(2) >0 && y(2) < 0), y = -y; end
    % Recall : p(0)=p(L)=0
     plot(x,y,'r-*') % solution par element finis
     title(['Mode ' num2str(m)]);
     xlabel('Position x/L');
     ylabel(' Normalized Amplitude ');
     legend('Analytical', 'FEM Quadratic ');
end