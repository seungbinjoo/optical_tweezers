%% Workspace initialisation

clear all;
close all;
clc;

%% Parameters

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

% number of nodes in each direction
l = 1;
p = 1;

%% Gaussian Beam

mu = [0, 0]; % Mean
sigma = [1, 0; 0, 1] * 2.6e-8; % Covariance matrix

% Generate grid
grid_lim = 3 * 1e-4;
grid_num_ticks = 20;
[x, y] = meshgrid(-grid_lim:(grid_lim/grid_num_ticks):grid_lim, -grid_lim:(grid_lim/grid_num_ticks):grid_lim);

% Combine x and y into a single matrix
xy = [x(:) y(:)];

% Calculate Gaussian function
z = mvnpdf(xy, mu, sigma);

% Reshape the result into a 2x2 matrix
relative_intensity = reshape(z, size(x));

% Normalize relative_intensity to [0, 1]
relative_intensity = relative_intensity / max(relative_intensity(:));

% Plot the Gaussian
z_shift = 5 * 1e-4;
gain = 15; % extra factor multipled to gaussian so the value heights are visible in plot

surf(x, y, (relative_intensity * 1e-5 * gain) - z_shift,'FaceAlpha', 1);
title('Optical Trap Simulation');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
axis equal
hold on

%% Plot rays before lens

% Define the endpoints of the vertical line
ray_before_length = 2 * 1e-4;

% Initialise vectors to store ray starting point and end point (xyz position values)
ray_before_start_point = zeros(numel(relative_intensity),3);
ray_before_end_point = zeros(numel(relative_intensity),3);

counter = 1;
for i = 1:size(x,1)
    for j = 1:size(y,1)
        ray_before_start_point(counter,1) = i*(grid_lim/grid_num_ticks) - grid_lim - (grid_lim/(2*grid_num_ticks));
        ray_before_start_point(counter,2) = j*(grid_lim/grid_num_ticks) - grid_lim - (grid_lim/(2*grid_num_ticks));
        ray_before_start_point(counter,3) = 0 - z_shift;
        
        ray_before_end_point(counter,1) = i*(grid_lim/grid_num_ticks) - grid_lim - (grid_lim/(2*grid_num_ticks));
        ray_before_end_point(counter,2) = j*(grid_lim/grid_num_ticks) - grid_lim - (grid_lim/(2*grid_num_ticks));
        ray_before_end_point(counter,3) = ray_before_length - z_shift;
        
        ray_before_x = [ray_before_start_point(counter,1) ray_before_end_point(counter,1)];
        ray_before_y = [ray_before_start_point(counter,2) ray_before_end_point(counter,2)];
        ray_before_z = [ray_before_start_point(counter,3) ray_before_end_point(counter,3)];
        
        % Plotting the 3D figure
        plot3(ray_before_x, ray_before_y, ray_before_z, 'r', 'LineWidth', 0.1); % Plotting the line
        grid on;
        hold on
        
        counter = counter + 1;
    end
end

% Set axes limits
xlim([-grid_lim, grid_lim]);
ylim([-grid_lim, grid_lim]);
zlim([-z_shift, 1 * 1e-4]); % Z AXES LIMITS!

%% Iris

% Define parameters for the annular ring
inner_radius = L;
outer_radius = 1;
height = ray_before_length;
num_points = 100;

% Generate points for the inner and outer circles
theta = linspace(0, 2*pi, num_points);
x_inner = inner_radius * cos(theta);
y_inner = inner_radius * sin(theta);
z_inner = ones(size(theta)) * height - z_shift; % Z-coordinates for inner circle

x_outer = outer_radius * cos(theta);
y_outer = outer_radius * sin(theta);
z_outer = ones(size(theta)) * height - z_shift; % Z-coordinates for outer circle

x_ring = [x_inner' x_outer'];
y_ring = [y_inner' y_outer'];
z_ring = [z_inner' z_inner'];

hold on
surf(x_ring, y_ring, z_ring,'FaceAlpha', 1)

% What are the intensity values from the Gaussian beam associated with each
% ray after the iris?
relative_intensity = reshape(relative_intensity,[numel(relative_intensity) 1]);

% eliminate rays that get blocked by the iris opening
counter = 1;

for i = 1:size(x,1)
    for j = 1:size(y,1)
        if sqrt(ray_before_end_point(counter,1)^2 + ray_before_end_point(counter,2)^2) > inner_radius
            ray_before_end_point(counter,1) = 0;
            ray_before_end_point(counter,2) = 0;
            ray_before_end_point(counter,3) = 0;
            relative_intensity(counter) = 0;
        end
        counter = counter + 1;
    end
end

% extend rays that haven't been eliminated
ray_after_iris_start_point = nonzeros(ray_before_end_point);
ray_after_iris_start_point = reshape(ray_after_iris_start_point, size(ray_after_iris_start_point,1)/3, 3);

ray_extend_after_iris = 2 * 1e-4;
ray_after_iris_end_point = ray_after_iris_start_point;
ray_after_iris_end_point(:,3) = ray_after_iris_start_point(:,3) + ray_extend_after_iris;

% What are the intensity values from the Gaussian beam associated with each
% ray after the iris?
relative_intensity_after_iris = nonzeros(relative_intensity);

% Plot new rays after iris
for i = 1:size(ray_after_iris_start_point,1)
    ray_after_iris_x = [ray_after_iris_start_point(i,1) ray_after_iris_end_point(i,1)];
    ray_after_iris_y = [ray_after_iris_start_point(i,2) ray_after_iris_end_point(i,2)];
    ray_after_iris_z = [ray_after_iris_start_point(i,3) ray_after_iris_end_point(i,3)];

    % Plotting the 3D figure
    plot3(ray_after_iris_x, ray_after_iris_y, ray_after_iris_z, 'r', 'LineWidth', 0.1); % Plotting the line
    grid on;
    hold on
end

%% Lens
NA = 1.20; % numerical aperture

ray_start_point = ray_after_iris_end_point;
ray_end_point = [zeros(size(ray_after_iris_end_point,1),2) (ray_before_length + height + f - z_shift)*ones(size(ray_after_iris_end_point,1),1)];

% Plot converging rays
for i = 1:size(ray_start_point,1)
    ray_x = [ray_start_point(i,1) ray_end_point(i,1)];
    ray_y = [ray_start_point(i,2) ray_end_point(i,2)];
    ray_z = [ray_start_point(i,3) ray_end_point(i,3)];

    % Plotting the 3D figure
    plot3(ray_x, ray_y, ray_z, 'r', 'LineWidth', 0.1); % Plotting the line
    grid on;
    hold on
end

%% Particle
[X,Y,Z] = sphere;

par_radius = 1e-6;
X2 = X * par_radius;
Y2 = Y * par_radius;
Z2 = Z * par_radius;

sphere_center = [0, 0, ray_before_length + height + f - z_shift];

surf(X2+sphere_center(1),Y2+sphere_center(2),Z2+sphere_center(3),'FaceAlpha', 0.5)

% Set transparency (alpha value) for the surface
alpha(0.5); % Adjust the value from 0 (fully transparent) to 1 (fully opaque)

hold on

%% Optical force calculation --> use OTGO

% PARAMETERS (AGAIN)
% Medium
nm = 1.33; % Medium refractive index --> water

% Spherical particle
R = 1e-6; % Particle radius [m] --> 1um radius particle
np = 1.59; % Particle refractive index --> polystyrene

% Focusing
f = 10e-6; % Focal length [m] --> ?
NA = 1.20; % numerical aperture --> water immersion objective
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m] --> ?
Ey0 = 1i*1e+4; % y electric field [V/m] --> ?
w0 = 20e-6; % Beam waist [m] --> ?
Nphi = 40; % Azimuthal divisions --> ?
Nr = 40; % Radial divisions --> ?
power = 10e-3; % Power [W] --> 10mW per microparticle for stable trapping

% Positions where to calculate the forces
x = [-1.6:0.32:1.6]*R;
y = [-1.6:0.32:1.6]*R;
z = [-2.4:0.32:2.4]*R;

% Initialization
% Trapping beam
blg = BeamLG(l,p,Ex0,Ey0,w0,R,Nphi,Nr);
blg = blg.normalize(power); % Set the power

% Calculates set of rays corresponding to focused optical beam
r = Ray.beam2focused(blg,f);

% PREALLOCATION OF VARIABLES
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

% For quiver3() function arguments
[X,Y,Z] = meshgrid(x,y,z);

% Transpose the 3D matrix (idk exactly why but this bug took me ages to solve)
X = permute(X,[2 1 3]);
Y = permute(Y,[2 1 3]);
Z = permute(Z,[2 1 3]);

% Plot 3D force vector field
quiver3(X,Y,Z,f_x,f_y,f_z,'b','LineWidth',2)

%% Plot 2D force field in transverse plane
% Extract 2D plane values where z is held constant
f_x_2D_tran = f_x(:,:,round(numel(z)/2));
f_y_2D_tran = f_y(:,:,round(numel(z)/2));

q_x_2D_tran = q_x(:,:,round(numel(z)/2));
q_y_2D_tran = q_y(:,:,round(numel(z)/2));

% Create 2D meshgrid
[X_tran, Y_tran] = meshgrid(x,y);

hold off; % release the hold
figure;

% Plot 2D force vector field
% Transpose the 2D matrices (idk exactly why but this bug took me ages to solve)
quiver(X_tran',Y_tran',f_x_2D_tran,f_y_2D_tran,'LineWidth',2)
title('2D Force Vector Field in Transverse Plane');
xlabel('X (m)');
ylabel('Y (m)');
xlim([-1.6 1.6]*R)
ylim([-1.6 1.6]*R)
axis equal
hold on

% Plot red circle
theta_circ = linspace(0, 2*pi, 100);  % Angle range for the circle
x_circle = R * cos(theta_circ);       % x-coordinates of the circle points
y_circle = R * sin(theta_circ);       % y-coordinates of the circle points
plot(x_circle, y_circle, 'LineWidth', 2);  % Plot the circle
xlim([-1.6 1.6]*R)
ylim([-1.6 1.6]*R)

%% Plot 2D force field in longitudinal plane
% Extract 2D plane values where y is held constant
f_x_2D_long = squeeze(f_x(:,round(numel(y)/2),:));
f_z_2D_long = squeeze(f_z(:,round(numel(y)/2),:));

q_x_2D_long = squeeze(q_x(:,round(numel(y)/2),:));
q_z_2D_long = squeeze(q_z(:,round(numel(y)/2),:));

% Create 2D meshgrid
[X_long, Z_long] = meshgrid(x,z);

hold off; % release the hold
figure;

% Plot 2D force vector field
% Transpose the 2D matrices (idk exactly why but this bug took me ages to solve)
quiver(X_long',Z_long',f_x_2D_long,f_z_2D_long,'LineWidth',2)
title('2D Force Vector Field in Longitudinal Plane');
xlabel('X (m)');
ylabel('Z (m)');
axis equal
xlim([-1.6 1.6]*R)
ylim([-2.4 2.4]*R)
hold on

% Plot red circle
theta_circ = linspace(0, 2*pi, 100);  % Angle range for the circle
x_circle = R * cos(theta_circ);       % x-coordinates of the circle points
y_circle = R * sin(theta_circ);       % y-coordinates of the circle points
plot(x_circle, y_circle, 'LineWidth', 2);  % Plot the circle