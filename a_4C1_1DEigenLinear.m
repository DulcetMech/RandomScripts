clear all; close all; clc;

% Problem data
Lx = 1; % domain length
c0 = 342; %speed of sound m/s

% step 1: Mesh
ne = input('Number of linear elements : '); % number of elements
nnt = ne + 1; % Total number of nodes
h = Lx/ne; % Length of the elements
x = [0:h:Lx];   % Coordinates table

%step 2 : Compute Elementary matrices
Ke = c0^2 * [1,-1;-1,1]*1/h;
Me = [2,1;1,2]*h/6; 

%Step 3 Assembling
I = eye(2,2);
K = zeros(nnt, nnt);
M = zeros(nnt, nnt);
for ie = 1:ne
    L = zeros(2, nnt);
    L(:, ie:ie + 1) = I;
    % Location matrixfor element ie
    K = K + L'*Ke*L;
    M = M + L'*Me*L;
end



% BC p(1)=p(nnt) = 0
K = K(2:nnt-1, 2:nnt-1);
M = M(2:nnt-1, 2:nnt-1);

ndof = nnt-2;

%compute eigen values and vectors

[V,D]=eig(K,M); D= sqrt(D);
Nor=V'*M*V; V=V*sqrt(inv(Nor)); % Normalization of the eigenvectors

for idof = 1: ndof
 w(idof) = D(idof,idof);
end

[w,ind]=sort(w); % Sort the eigenvalues in ascending order


% Step 6: Comparison with the exact solution for a selected number of modes

nm = min(4, ndof); % Here we select the first modes up to the fourth
theo = c0*[1:ndof]*pi / Lx;
Err = abs(theo-w)./theo;
fprintf('\n ************Results for %d lin elements ************', ne);
fprintf('\n ====================================================');
fprintf('\n Mode FEM (Hz) Exact (Hz) Error(%%)');
fprintf('\n --------------------------------');
for m=1:nm
 fprintf('\n %d %7.2f %7.2f %3.1f', m, w(m)/2/pi, theo(m)/2/pi, Err(m)*100);
end
fprintf('\n ================================================\n');
fprintf('\n');
%
% comparison of the mode shapes
%
xt=[0:0.025:1] * Lx; % Mesh to visualize the mode shapes
for m=1:nm
 figure;
 yt= sqrt(2)*sin(m*pi*xt);
 plot(xt,yt,'LineWidth',2,'Color',[0 0 0])
 hold on, % Theoretical solution for

 y=[0,real(V(:,ind(m))'),0]; if (yt(2) >0 && y(2) < 0), y = -y; end;
% Recall : p(0)=p(L)=0
 plot(x,y,'-*','LineWidth',2,'Color',[0.5 0.5 0.5]) % FE solution
 title(['Mode' num2str(m)]);
 xlabel('Position x/L');
 ylabel(' Normalized Amplitude ');
 legend('Analytical', 'FEM Linear');
end


