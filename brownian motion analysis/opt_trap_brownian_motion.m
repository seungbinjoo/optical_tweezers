%% Workspace initialisation

clear all;
close all;
clc;

%% 1. Brownian dynamics simulation
%% 1.1. Parameters for the optical trap
% Medium
nm = 1.33; % Medium refractive index --> water

% Spherical particle
R = 1e-6; % Particle radius [m] --> 1um radius particle
np = 1.54; % Particle refractive index --> polystyrene

% Focusing
f = 2.525e-3; % Focal length [m]
NA = 1.1; % numerical aperture --> water immersion objective
L = f*NA/nm; % Iris aperture radius [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m] --> linearly polarised
Ey0 = 0; % y electric field [V/m] --> linearly polarised
w0 = L; % Beam waist [m] = beam radius at lens = overfilling so 1/e^2 matches aperture edge
Nphi = 40; % Azimuthal divisions 
Nr = 40; % Radial division
power = 10e-3; % Power [W] --> 10 mW

% Initialization
% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to focused optical beam
r = Ray.beam2focused(bg,f);

%% 1.2. Parameters for brownian motion simulation

% Brownian motion
kB = PhysConst.kB; % Boltzmann constant, 1.3806e-23 J/K
T = 293; % Temperature [K] --> 20 degrees Celsius ()
eta = 0.001; % Water viscosity [Pa*s]
time_step = 1e-3; % timestep [s]
iterations = 500; % number of steps

% Diffusion Constant of the Spherical Particle
gamma = (6*pi*eta*R);
D = kB*T / gamma;

% Define trap stiffnesses in each direction
kappa_x = 5.83e-12 / 1e-6; % 5.83 pN/um
kappa_y = 5.83e-12 / 1e-6; % 5.83 pN/um
kappa_z = 0.83e-12 / 1e-6; % 0.83 pN/um

% Number of trials
trials = 4;

% Initial position of the center of the microparticle for each trial
x_init = [1*R -1*R 1*R -1*R];
y_init = [1*R 1*R -1*R -1*R];
z_init = [0.5*R 0.5*R 0.5*R 0.5*R];


%% 1.3. Simulation loop (Euler's method and brownian dynamics)
% Preallocation of matrix that will store particle's center positions (x,y,z)
position = zeros(3*trials,iterations);

% Dynamics is based on page 202 equation 7.24 in optical tweezers textbook
for k = 1:trials

    % Initial position for the trial
    x = x_init(k);
    y = y_init(k);
    z = z_init(k);

    for i = 1:iterations
        
        position((k-1)*3+1,i) = x;
        position((k-1)*3+2,i) = y;
        position((k-1)*3+3,i) = z;

        % To track progress of the simulation
        disp(['current time = ' num2str(i*time_step) ' s / ' num2str(iterations*time_step) ' s'])
    
        % % Spherical particle
        % bead = ParticleSpherical(Point(x,y,z),R,nm,np);
        % 
        % forces = bead.force(r);
        % force = Vector(x,y,z, ...
        %     sum(forces.Vx(isfinite(forces.Vx))), ...
        %     sum(forces.Vy(isfinite(forces.Vy))), ...
        %     sum(forces.Vz(isfinite(forces.Vz))) ...
        %     );
        
        % Particle position update - optical force contribution
        x = x - (time_step*kappa_x*x)/gamma;
        y = y - (time_step*kappa_y*y)/gamma;
        z = z - (time_step*kappa_z*z)/gamma;
        
        
        % Particle position update - brownian motion contribution
        x = x + sqrt(2*time_step*D)*randn();
        y = y + sqrt(2*time_step*D)*randn();
        z = z + sqrt(2*time_step*D)*randn();
    
    end
end

%% 1.4. Plot simulation results - first half

x_plot_1 = position(1,1:iterations/2);
y_plot_1 = position(2,1:iterations/2);
z_plot_1 = position(3,1:iterations/2);

x_plot_2 = position(4,1:iterations/2);
y_plot_2 = position(5,1:iterations/2);
z_plot_2 = position(6,1:iterations/2);

x_plot_3 = position(7,1:iterations/2);
y_plot_3 = position(8,1:iterations/2);
z_plot_3 = position(9,1:iterations/2);

x_plot_4 = position(10,1:iterations/2);
y_plot_4 = position(11,1:iterations/2);
z_plot_4 = position(12,1:iterations/2);

tiledlayout(2,trials,"TileSpacing","tight","Padding",'compact')
nexttile([1 trials])

plot3(x_plot_1,y_plot_1,z_plot_1,'LineWidth',1)
grid
axis equal
title('Brownian motion in optical trap, trial 1-4, first 250 ms','FontSize',14)
xlabel('x (m)','FontSize',14)
ylabel('y (m)','FontSize',14)
zlabel('z (m)','FontSize',14)
set(gca, 'FontSize', 12)
hold on
plot3(x_plot_2,y_plot_2,z_plot_2,'LineWidth',1)
hold on
plot3(x_plot_3,y_plot_3,z_plot_3,'LineWidth',1)
hold on
plot3(x_plot_4,y_plot_4,z_plot_4,'LineWidth',1)
hold off

legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4','FontSize',14)
legend('Location', 'southeast');

% Plot sphere representing the volume in which the particle center needs to
% be inside in order for the particle to be succesfully inserted
% abl_hole_diameter = 4e-6;
% R_center = (abl_hole_diameter - 2*R)/2; % radius of the sphere in which the microparticle center needs to be in order for successful insertion
% 
% [x_sphere, y_sphere, z_sphere] = sphere; % Generate sphere coordinates
% x_sphere_scaled = R_center * x_sphere + 0; % Scale and translate sphere coordinates
% y_sphere_scaled = R_center * y_sphere + 0;
% z_sphere_scaled = R_center * z_sphere + 0.3e-6; % trap center is 0.3 um in the positive z direction
% surface(x_sphere_scaled, y_sphere_scaled, z_sphere_scaled,'FaceColor', [0.9290 0.6940 0.1250],'FaceAlpha', 0.1,'EdgeColor',[0.8 0.8 0.8])
% 
% hold off

%% 1.5. Plot simulation results - second half
x_plot_1 = position(1,iterations/2+1:iterations);
y_plot_1 = position(2,iterations/2+1:iterations);
z_plot_1 = position(3,iterations/2+1:iterations);

x_plot_2 = position(4,iterations/2+1:iterations);
y_plot_2 = position(5,iterations/2+1:iterations);
z_plot_2 = position(6,iterations/2+1:iterations);

x_plot_3 = position(7,iterations/2+1:iterations);
y_plot_3 = position(8,iterations/2+1:iterations);
z_plot_3 = position(9,iterations/2+1:iterations);

x_plot_4 = position(10,iterations/2+1:iterations);
y_plot_4 = position(11,iterations/2+1:iterations);
z_plot_4 = position(12,iterations/2+1:iterations);

nexttile
plot3(x_plot_1,y_plot_1,z_plot_1,'LineWidth',1,'Color',[0 0.4470 0.7410])
grid
axis equal
title('Trial 1, last 250 ms','FontSize',14)
xlabel('x (m)','FontSize',14)
ylabel('y (m)','FontSize',14)
zlabel('z (m)','FontSize',14)
set(gca, 'FontSize', 12)

nexttile
plot3(x_plot_2,y_plot_2,z_plot_2,'LineWidth',1,'Color',[0.8500 0.3250 0.0980])
grid
axis equal
title('Trial 2, last 250 ms','FontSize',14)
xlabel('x (m)','FontSize',14)
ylabel('y (m)','FontSize',14)
zlabel('z (m)','FontSize',14)
set(gca, 'FontSize', 12)

nexttile
plot3(x_plot_3,y_plot_3,z_plot_3,'LineWidth',1,'Color',[0.9290 0.6940 0.1250])
grid
axis equal
title('Trial 3, last 250 ms','FontSize',14)
xlabel('x (m)','FontSize',14)
ylabel('y (m)','FontSize',14)
zlabel('z (m)','FontSize',14)
set(gca, 'FontSize', 12)

nexttile
plot3(x_plot_4,y_plot_4,z_plot_4,'LineWidth',1,'Color',[0.4940 0.1840 0.5560])
grid
axis equal
title('Trial 4, last 250 ms','FontSize',14)
xlabel('x (m)','FontSize',14)
ylabel('y (m)','FontSize',14)
zlabel('z (m)','FontSize',14)
set(gca, 'FontSize', 12)

% Plot sphere representing the volume in which the particle center needs to
% be inside in order for the particle to be succesfully inserted
% abl_hole_diameter = 4e-6;
% R_center = (abl_hole_diameter - 2*R)/2; % radius of the sphere in which the microparticle center needs to be in order for successful insertion
% 
% [x_sphere, y_sphere, z_sphere] = sphere; % Generate sphere coordinates
% x_sphere_scaled = R_center * x_sphere + 0; % Scale and translate sphere coordinates
% y_sphere_scaled = R_center * y_sphere + 0;
% z_sphere_scaled = R_center * z_sphere + 0.3e-6; % trap center is 0.3 um in the positive z direction
% surface(x_sphere_scaled, y_sphere_scaled, z_sphere_scaled,'FaceColor', [0.9290 0.6940 0.1250],'FaceAlpha', 0.1,'EdgeColor',[0.8 0.8 0.8])
% 
% hold off
