% *********************************************************************** %
% Lattice Boltzmann Method (LBM)                                          %
%                                                                         %
% Email: zhe.li@ec-nantes.fr                                              %
% Version: 0.0                                                            %
% Lattice: D2Q9                                                           %
%                         6---2---5                                       %
%                         | \ | / |                                       %
%                         3---0---1                                       %
%                         | / | \ |                                       %
%                         7---4---8                                       %
%                                                                         %
% *********************************************************************** %
%%
clear variables
close all
clc
format long
%%
%
% ======================================================================= %
%                               Constants                                 %
% ======================================================================= %
cs_lat = 1.0/(3^0.5);
cs2 = cs_lat^2.0;
cs4 = cs2^2.0;
inv1cs2 = 1.0/cs2;
inv2cs2 = 1.0/(2.0*cs2);
inv1cs4 = 1.0/cs4;
inv2cs4 = 1.0/(2.0*cs4);
ex = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0]';
ey = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]';
w = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, ...
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0];
nbDist = 9;


% ======================================================================= %
%                       Computational parameters                          %
% ======================================================================= %
Nx = 3; 
Ny = 21; 
nbSteps = 1000; % number of time steps

L = 10; % periodic physical domain horizontal length [m]
H = 3;  % periodic physical domain height [m]
time = 5; % simulation time [s]

Cx  = L/(Nx-1); % horizontal length scaling factor [m]
Ct  = time/nbSteps; % time scaling factor [s]

Cu  = Cx/Ct; % velocity scaling factor [m/s]
Cnu = Cx^2/Ct; % kinematic viscosity scaling factor [m^2/s]
Cg = Cx/Ct^2;

nu = 100; % kinematic viscosity [m^2/s]
nu_LBM = nu/Cnu; 

tau = nu_LBM/cs2 + .5; % relaxation time

g   = [1; 0]./Cg; % constant horizontal body force [m/s^s]

% check values of Ma, Re and tau
H_LBM    = H/Cx;
Umax_LBM = g(1)*H^2/8/nu_LBM; % maximum velocity at y = .5*H m LBM
Ma       = Umax_LBM/cs_lat; % maximum Mach number

Umean_LBM = g(1)*H_LBM/12/nu_LBM; % mean velocity LBM
Re        = Umean_LBM*H_LBM/nu_LBM; % Reynolds number

if Ma > 0.6
    error('Ma < 0.6, compressible flow !')

end

if abs(tau - 0.5) < 0.05
    error('Instability due to too low relaxation time for BGK !')

end

if Re > 2100
    error('Re > 2100, turbulent flow')

end

% ======================================================================= %
%                           Initialization                                %
% ======================================================================= %
fdist = zeros(nbDist, Nx, Ny); % distribution function matrix
fcoll = zeros(nbDist, Nx, Ny); % post-collision distribution function matrix
rho   = ones(Nx, Ny);
vx    = zeros(Nx, Ny);
vy    = zeros(Nx, Ny);

Jx0   = rho.*vx;
Jy0   = rho.*vy;
rho0  = rho;

for j = 1:Ny
    for i = 1:Nx
        u = [vx(i,j); vy(i,j)];
        for k = 1:nbDist
            csi = [ex(k);ey(k)];
            fdist(k,i,j) = rho(i,j)*w(k)*( 1 + inv1cs2*csi'*u + inv2cs4*(csi'*u)^2 - inv2cs2*u'*u);
        end
    end
end

% ======================================================================= %
%                             Iterations                                  %
% ======================================================================= %
for iStep = 1:nbSteps
    % ---------------------------------
    % (1) Collision
    for j = 1:Ny
        for i = 1:Nx
            u = [vx(i,j);vy(i,j)];
            for k = 1:nbDist
                csi = [ex(k);ey(k)];
                Fa  = rho(i,j)*w(k)*( inv1cs2*(csi - u) + inv1cs4*(csi'*u)*csi )'*g;
                feq = rho(i,j)*w(k)*( 1 + inv1cs2*csi'*u + inv2cs4*(csi'*u)^2 - inv2cs2*u'*u);
                fcoll(k,i,j) = fdist(k,i,j) - (fdist(k,i,j) - feq)/tau + Fa;
            end
        end
    end
    
    % ---------------------------------
    % (2) Streaming
    for j = 1:Ny
        for i = 1:Nx
            for k = 1:nbDist
                csi = [ex(k);ey(k)];
                ia = i + csi(1);
                ja = j + csi(2);
                if ja > 0 && ja <= Ny

                    if ia > 0 && ia <= Nx
                        fdist(k,ia,ja) = fcoll(k,i,j);
                    end

                    % periodic boundary condition
                    if ia < 1
                        ia = Nx;
                        fdist(k,ia,ja) = fcoll(k,i,j);
                    end
                    if ia > Nx
                        ia = 1;
                        fdist(k,ia,ja) = fcoll(k,i,j);
                    end

                end
            end
        end
    end
    
    % ---------------------------------
    % (3) Boundary condition
    % (3.3) Bottom
    j = 1;
    for i = 1:Nx
        fdist(3,i,j) = fcoll(5,i,j);
        fdist(6,i,j) = fcoll(8,i,j);
        fdist(7,i,j) = fcoll(9,i,j);
    end

    % (3.4) Top
    j = Ny;
    for i = 1:Nx
        fdist(5,i,j) = fcoll(3,i,j);
        fdist(8,i,j) = fcoll(6,i,j);
        fdist(9,i,j) = fcoll(7,i,j);
    end
    
    % ---------------------------------
    % (4) Macroscopic variables
    residu = 0.0;
    for j = 1:Ny
        for i = 1:Nx
            rho(i,j) = sum(fdist(:,i,j));
            Jx(i,j)  = sum(fdist(:,i,j).*ex);
            Jy(i,j)  = sum(fdist(:,i,j).*ey);

            Vx(i,j)  = Jx(i,j)/rho(i,j);
            Vy(i,j)  = Jy(i,j)/rho(i,j);

            residu = residu + abs(rho(i,j) - rho0(i,j))/rho0(i,j);
        end
    end
    fprintf('\n Mass conservation relative residual res = %g \n', residu)
end

%%

figure
hold on
grid on
ylim([0 Ny])
plot(Vx(2,:),[1:1:Ny],'bo-')
hold off
