%% Workspace initialisation

clear all;
close all;
clc;

%% Simulation parameters - beam waist = 0.3 x (iris aperture)
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
w0 = 0.3*L; % Beam waist [m] = beam radius at lens
Nphi = 40; % Azimuthal divisions 
Nr = 40; % Radial division
power = 5e-3; % Power [W] --> 5 mW

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
r = Ray.beam2focused(bg,f);

%% Optical trap simulation - transverse plane
% PREALLOCATION OF VARIABLES
% Force components
f_x = zeros(numel(x), numel(y));
f_y = zeros(numel(x), numel(y));
f_z = zeros(numel(x), numel(y));

% Force magnitude
f_tran = zeros(numel(x), numel(z));

% Trapping efficiencies
q_x = zeros(numel(x), numel(y));
q_y = zeros(numel(x), numel(y));
q_z = zeros(numel(x), numel(y));

% Trapping efficiency based on force magnitude
q_tran = zeros(numel(x), numel(y));

% TRAPPING EFFICIENCY CALCULATION
for i = 1:1:length(x)
    for j = 1:1:length(y)
    
            % Display update message
            disp(['Trapping efficiency calculation for transverse plane ' int2str((i-1)*length(z)+j) '/' int2str(numel(f_x))])

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

            % Trapping efficiencies in each direction
            q_x(i,j) = f_x(i,j)/power*(PhysConst.c0/nm);
            q_y(i,j) = f_y(i,j)/power*(PhysConst.c0/nm);
            q_z(i,j) = f_z(i,j)/power*(PhysConst.c0/nm);

            % Trapping efficiency in transverse plane
            f_tran(i,j) = sqrt(f_x(i,j).^2 + f_y(i,j).^2);
            q_tran(i,j) = f_tran(i,j)/power*(PhysConst.c0/nm);
    
    end
end

%% Plot trapping efficiency surface in transverse plane
% Create 2D meshgrid
[X_tran, Y_tran] = meshgrid(x,y);

tiledlayout(2,2,"TileSpacing","compact")
nexttile

% Plot 2D trapping efficiency field
% Transpose the 2D matrices (idk exactly why but this bug took me ages to solve)
surf(Y_tran',X_tran',abs(q_tran))
title('Trapping Efficiency in Transverse Plane (w_0 = 0.3L)','FontSize',14);
xlabel('Y (m)','FontSize',16);
ylabel('X (m)','FontSize',16);
zlabel('Trapping efficiency')
set(gca, 'FontSize', 16)
view(2)
colorbar

%% Optical trap simulation - longitudinal plane
% PREALLOCATION OF VARIABLES
% Force components
f_x = zeros(numel(x), numel(z));
f_y = zeros(numel(x), numel(z));
f_z = zeros(numel(x), numel(z));

% Force magnitude
f_long = zeros(numel(x), numel(z));

% Trapping efficiencies
q_x = zeros(numel(x), numel(z));
q_y = zeros(numel(x), numel(z));
q_z = zeros(numel(x), numel(z));

% Trapping efficiency based on force magnitude
q_long = zeros(numel(x), numel(z));

% TRAPPING EFFICIENCY CALCULATION
for i = 1:1:length(x)
    for j = 1:1:length(z)
    
            % Display update message
            disp(['Trapping efficiency calculation for longitudinal plane ' int2str((i-1)*length(z)+j) '/' int2str(numel(f_x))])

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

            % Trapping efficiencies in each direction
            q_x(i,j) = f_x(i,j)/power*(PhysConst.c0/nm);
            q_y(i,j) = f_y(i,j)/power*(PhysConst.c0/nm);
            q_z(i,j) = f_z(i,j)/power*(PhysConst.c0/nm);

            % Trapping efficiency in longitudinal plane
            f_long(i,j) = sqrt(f_x(i,j).^2 + f_z(i,j).^2);
            q_long(i,j) = f_long(i,j)/power*(PhysConst.c0/nm);
    
    end
end

%% Plot trapping efficiency surface in longitudinal plane
% Create 2D meshgrid
[X_long, Z_long] = meshgrid(x,z);

nexttile

% Plot trapping effiency surface
% Transpose the 2D matrices (idk exactly why but this bug took me ages to solve)
surf(Z_long',X_long',abs(q_long))
title('Trapping Efficiency in Longitudinal Plane (w_0 = 0.3L)','FontSize',14);
xlabel('Z (m)','FontSize',16);
ylabel('X (m)','FontSize',16);
zlabel('Trapping efficiency')
set(gca, 'FontSize', 16)
view(2)
colorbar

%% Simulation parameters - beam waist = 1.0 x (iris aperture)
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
power = 5e-3; % Power [W] --> 5 mW

% Positions where to calculate the forces
x = [-1.6:increment:1.6]*R;
y = [-1.6:increment:1.6]*R;
z = [-1.6:increment:1.6]*R;

% Initialization
% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to focused optical beam
r = Ray.beam2focused(bg,f);

%% Optical trap simulation - transverse plane
% PREALLOCATION OF VARIABLES
% Force components
f_x = zeros(numel(x), numel(y));
f_y = zeros(numel(x), numel(y));
f_z = zeros(numel(x), numel(y));

% Force magnitude
f_tran = zeros(numel(x), numel(z));

% Trapping efficiencies
q_x = zeros(numel(x), numel(y));
q_y = zeros(numel(x), numel(y));
q_z = zeros(numel(x), numel(y));

% Trapping efficiency based on force magnitude
q_tran = zeros(numel(x), numel(y));

% TRAPPING EFFICIENCY CALCULATION
for i = 1:1:length(x)
    for j = 1:1:length(y)
    
            % Display update message
            disp(['Trapping efficiency calculation for transverse plane ' int2str((i-1)*length(z)+j) '/' int2str(numel(f_x))])

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

            % Trapping efficiencies in each direction
            q_x(i,j) = f_x(i,j)/power*(PhysConst.c0/nm);
            q_y(i,j) = f_y(i,j)/power*(PhysConst.c0/nm);
            q_z(i,j) = f_z(i,j)/power*(PhysConst.c0/nm);

            % Trapping efficiency in transverse plane
            f_tran(i,j) = sqrt(f_x(i,j).^2 + f_y(i,j).^2);
            q_tran(i,j) = f_tran(i,j)/power*(PhysConst.c0/nm);
    
    end
end

%% Plot trapping efficiency surface in transverse plane
% Create 2D meshgrid
[X_tran, Y_tran] = meshgrid(x,y);

nexttile

% Plot 2D trapping efficiency field
% Transpose the 2D matrices (idk exactly why but this bug took me ages to solve)
surf(Y_tran',X_tran',abs(q_tran))
title('Trapping Efficiency in Transverse Plane (w_0 = L)', 'FontSize', 14);
xlabel('Y (m)', 'FontSize', 16);
ylabel('X (m)', 'FontSize', 16);
zlabel('Trapping efficiency')
set(gca, 'FontSize', 16)
view(2)
colorbar

%% Optical trap simulation - longitudinal plane
% PREALLOCATION OF VARIABLES
% Force components
f_x = zeros(numel(x), numel(z));
f_y = zeros(numel(x), numel(z));
f_z = zeros(numel(x), numel(z));

% Force magnitude
f_long = zeros(numel(x), numel(z));

% Trapping efficiencies
q_x = zeros(numel(x), numel(z));
q_y = zeros(numel(x), numel(z));
q_z = zeros(numel(x), numel(z));

% Trapping efficiency based on force magnitude
q_long = zeros(numel(x), numel(z));

% TRAPPING EFFICIENCY CALCULATION
for i = 1:1:length(x)
    for j = 1:1:length(z)
    
            % Display update message
            disp(['Trapping efficiency calculation for longitudinal plane ' int2str((i-1)*length(z)+j) '/' int2str(numel(f_x))])

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

            % Trapping efficiencies in each direction
            q_x(i,j) = f_x(i,j)/power*(PhysConst.c0/nm);
            q_y(i,j) = f_y(i,j)/power*(PhysConst.c0/nm);
            q_z(i,j) = f_z(i,j)/power*(PhysConst.c0/nm);

            % Trapping efficiency in longitudinal plane
            f_long(i,j) = sqrt(f_x(i,j).^2 + f_z(i,j).^2);
            q_long(i,j) = f_long(i,j)/power*(PhysConst.c0/nm);
    
    end
end

%% Plot trapping efficiency surface in longitudinal plane
% Create 2D meshgrid
[X_long, Z_long] = meshgrid(x,z);

nexttile

% Plot trapping effiency surface
% Transpose the 2D matrices (idk exactly why but this bug took me ages to solve)
surf(Z_long',X_long',abs(q_long))
title('Trapping Efficiency in Longitudinal Plane (w_0 = L)', 'FontSize', 14);
xlabel('Z (m)', 'FontSize', 16);
ylabel('X (m)', 'FontSize', 16);
zlabel('Trapping efficiency')
set(gca, 'FontSize', 16)
view(2)
colorbar