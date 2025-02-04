%% Workspace initialisation

clear all;
close all;
clc;

%% Simulation parameters - beam waist = 2 x (iris aperture)
% Medium
nm = 1.33; % Medium refractive index --> water

% Spherical particle
R = 1e-6; % Particle radius [m] --> 1um radius particle
np = 1.59; % Particle refractive index --> polystyrene

% Focusing
f_long = 100e-6; % Focal length [m] --> ?
NA = 1.20; % numerical aperture --> water immersion objective
L = f_long*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m] --> ?
Ey0 = 1i*1e+4; % y electric field [V/m] --> ?
w0 = 20e-6; % Beam waist [m] --> ?
Nphi = 40; % Azimuthal divisions --> ?
Nr = 40; % Radial divisions --> ?
power = 10e-3; % Power [W] --> 10mW per microparticle for stable trapping

% Increment (i.e. size of grid ticks in 2D plots)
increment = 0.08;

% Positions where to calculate the forces
x = [-1.6:increment:1.6]*R;
y = [-1.6:increment:1.6]*R;
z = [-1.6:increment:1.6]*R;

% Initialization
% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to focused optical beam
r = Ray.beam2focused(bg,f_long);

%% Optical trap simulation - longitudinal plane
% PREALLOCATION OF VARIABLES
% Force components
f_x = zeros(numel(x), numel(z));
f_y = zeros(numel(x), numel(z));
f_z = zeros(numel(x), numel(z));

% Force magnitude
f_long = zeros(numel(x), numel(z));

% Trapping stiffnesses
kappa_x = zeros(numel(x), numel(z));
kappa_y = zeros(numel(x), numel(z));
kappa_z = zeros(numel(x), numel(z));

% Trapping stiffness based on force magnitude
kappa_long = zeros(numel(x), numel(z));

% TRAPPING STIFFNESS CALCULATION
for i = 1:1:length(x)
    for j = 1:1:length(z)
    
            % Display update message
            disp(['Trap stiffness calculation for longitudinal plane ' int2str((i-1)*length(z)+j) '/' int2str(numel(f_x))])

            % Spherical particle
            bead = ParticleSpherical(Point(x(i),0,z(j)),R,nm,np);

            % Calculate force
            forces = bead.force(r);
            force = Vector(bead.sp.c.X,bead.sp.c.Y,bead.sp.c.Z, ...
                sum(forces.Vx(isfinite(forces.Vx))), ...
                sum(forces.Vy(isfinite(forces.Vy))), ...
                sum(forces.Vz(isfinite(forces.Vz))) ...
                );

            % Trapping optical force in each direction
            f_x(i,j) = force.Vx;
            f_y(i,j) = force.Vy;
            f_z(i,j) = force.Vz;

            % Trapping stiffness in each direction
            kappa_x(i,j) = f_x(i,j)/power*(PhysConst.c0/nm);
            kappa_y(i,j) = f_y(i,j)/power*(PhysConst.c0/nm);
            kappa_z(i,j) = f_z(i,j)/power*(PhysConst.c0/nm);

            % Trapping stiffness in longitudinal plane
            f_long(i,j) = sqrt(f_x(i,j).^2 + f_z(i,j).^2);
            kappa_long(i,j) = f_long(i,j)/power*(PhysConst.c0/nm);
    
    end
end

%% Plot trapping stiffness surface in longitudinal plane
% Create 2D meshgrid
[X_long, Z_long] = meshgrid(x,z);

tiledlayout(2,2,"TileSpacing","compact")
nexttile

% Plot trapping stiffness surface
% Transpose the 2D matrices (idk exactly why but this bug took me ages to solve)
surf(Z_long',X_long',abs(kappa_long))
title('2D Trap stiffness Field in Longitudinal Plane (with overfilling)');
xlabel('Z (m)');
ylabel('X (m)');
zlabel('Trap stiffness')
view(2)

%% Optical trap simulation - transverse plane
% PREALLOCATION OF VARIABLES
% Force components
f_x = zeros(numel(x), numel(y));
f_y = zeros(numel(x), numel(y));
f_z = zeros(numel(x), numel(y));

% Force magnitude
f_tran = zeros(numel(x), numel(z));

% Trapping stiffness
kappa_x = zeros(numel(x), numel(y));
kappa_y = zeros(numel(x), numel(y));
kappa_z = zeros(numel(x), numel(y));

% Trapping stiffness based on force magnitude
kappa_tran = zeros(numel(x), numel(y));

% TRAPPING STIFFNESS CALCULATION
for i = 1:1:length(x)
    for j = 1:1:length(y)
    
            % Display update message
            disp(['Trap stiffness calculation for transverse plane ' int2str((i-1)*length(z)+j) '/' int2str(numel(f_x))])

            % Spherical particle
            bead = ParticleSpherical(Point(x(i),y(j),0),R,nm,np);

            % Calculate force
            forces = bead.force(r);
            force = Vector(bead.sp.c.X,bead.sp.c.Y,bead.sp.c.Z, ...
                sum(forces.Vx(isfinite(forces.Vx))), ...
                sum(forces.Vy(isfinite(forces.Vy))), ...
                sum(forces.Vz(isfinite(forces.Vz))) ...
                );

            % Trapping optical force in each direction
            f_x(i,j) = force.Vx;
            f_y(i,j) = force.Vy;
            f_z(i,j) = force.Vz;

            % Trapping stiffness in each direction
            kappa_x(i,j) = f_x(i,j)/power*(PhysConst.c0/nm);
            kappa_y(i,j) = f_y(i,j)/power*(PhysConst.c0/nm);
            kappa_z(i,j) = f_z(i,j)/power*(PhysConst.c0/nm);

            % Trapping stiffness in transverse plane
            f_tran(i,j) = sqrt(f_x(i,j).^2 + f_y(i,j).^2);
            kappa_tran(i,j) = f_tran(i,j)/power*(PhysConst.c0/nm);
    
    end
end

%% Plot trapping stiffness surface in transverse plane
% Create 2D meshgrid
[X_tran, Y_tran] = meshgrid(x,y);

nexttile

% Plot 2D trapping stiffness field
% Transpose the 2D matrices (idk exactly why but this bug took me ages to solve)
surf(Y_tran',X_tran',abs(kappa_tran))
title('2D Trap Stiffness Field in Transverse Plane (with overfilling)');
xlabel('Y (m)');
ylabel('X (m)');
zlabel('Trap stiffness')
view(2)