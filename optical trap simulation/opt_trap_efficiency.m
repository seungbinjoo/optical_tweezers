%% Workspace initialisation

clear all;
close all;
clc;

%% Define parameters for optical trapping

% Medium
nm = 1.33; % Medium refractive index --> water

% Spherical particle
R = 1e-6; % Particle radius [m] --> 1um radius particle
np = 1.59; % Particle refractive index --> polystyrene

% Focusing
f = 100e-6; % Focal length [m] --> ?
NA = 1.20; % numerical aperture --> water immersion objective
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m] --> ?
Ey0 = 1i*1e+4; % y electric field [V/m] --> ?
w0 = 100e-6; % Beam waist [m] --> ?
Nphi = 40; % Azimuthal divisions --> ?
Nr = 40; % Radial divisions --> ?
power = 10e-3; % Power [W] --> 10mW per microparticle for stable trapping

% Positions along axis where to calculate the trapping efficiencies
x = [-1.6:0.32:1.6]*R; % transversal axis
z = [-1.6:0.32:1.6]*R; % longitudinal axis

%% PART 1: TRANSVERSAL AXIS
%% Initialization of the incoming Gaussian beam
% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to focused optical beam
r = Ray.beam2focused(bg,f);

%% Preallocation of variables that will store trapping forces and efficiencies
% Force components
f_x = zeros(numel(x), numel(y), numel(z));
f_y = zeros(numel(x), numel(y), numel(z));
f_z = zeros(numel(x), numel(y), numel(z));

% Trapping efficiencies
q_x = zeros(numel(x), numel(y), numel(z));
q_y = zeros(numel(x), numel(y), numel(z));
q_z = zeros(numel(x), numel(y), numel(z));

% FORCE CALCULATION
for i = 1:1:length(x)
    for j = 1:1:length(y)
        for k = 1:1:length(z)
    
            % Display update message
            disp(['Force calculation ' int2str((i-1)*length(y)*length(z)+(j-1)*length(z)+k) '/' int2str(numel(f_x))])

            % Spherical particle
            bead = ParticleSpherical(Point(x(i),y(j),z(k)),R,nm,np);

            % Calculate force
            forces = bead.force(r);
            force = Vector(bead.sp.c.X,bead.sp.c.Y,bead.sp.c.Z, ...
                sum(forces.Vx(isfinite(forces.Vx))), ...
                sum(forces.Vy(isfinite(forces.Vy))), ...
                sum(forces.Vz(isfinite(forces.Vz))) ...
                );

            % Trapping optical force in each direction
            f_x(i,j,k) = force.Vx;
            f_y(i,j,k) = force.Vy;
            f_z(i,j,k) = force.Vz;

            % Trapping efficiencies in each direction
            q_x(i,j,k) = f_x(i,j,k)/power*(PhysConst.c0/nm);
            q_y(i,j,k) = f_y(i,j,k)/power*(PhysConst.c0/nm);
            q_z(i,j,k) = f_z(i,j,k)/power*(PhysConst.c0/nm);
    
        end
    end
end

