%% Workspace initialization
clear all;
close all;
clc;

%% 1. Generate E_beam (incoming beam)
% Parameters
lambda_beam = 975 * 10e-9;     % Wavelength = 975 nm
w0_beam = 17.4 * 1e-3;         % Beam diameter (beam waist) = 17.40 mm (from LCOS SLM specifications)
z_beam = 0;                    % Propagation distance (in meters)
shift = [0 0];                 % Translating the beam in xy space
resolution = 1024;

% SLM dimensions
x_slm = linspace((-17.4/2)*1e-3,(17.4/2)*1e-3,resolution);
y_slm = linspace((-17.4/2)*1e-3,(17.4/2)*1e-3,resolution);

% Create a grid
x_beam = linspace(-5*w0_beam, 5*w0_beam, resolution);
y_beam = linspace(-5*w0_beam, 5*w0_beam, resolution);

% Compute intensity of single Gaussian beam (incoming light from the expanded laser light)
E_beam = computeGaussian(lambda_beam, w0_beam, z_beam, shift, x_slm, y_slm);

% Plot the Gaussian beam profile
figure(1);
imagesc(x_slm, y_slm, abs(E_beam).^2);
axis xy;
colormap('default');
title('Incoming Gaussian Beam Profile E\_{beam}');
xlabel('X (m)');
ylabel('Y (m)');
colorbar;


%% 2. Generate E_target (target/desired distribution at the focal plane)
% Parameters
lambda_target = 975 * 1e-9;     % Wavelength = 975 nm
w0_target = 0.74 * 1e-6;        % Beam radius at focus = 0.37 um
z_target = 0;                   % Propagation distance (in meters)
n_particles = 1;                % Number of particles we want to tweeze = number of beams we want to generate

% Create a grid
FOV = 150e-6;
x_target = linspace(-FOV/2, FOV/2, resolution);
y_target = linspace(-FOV/2, FOV/2, resolution);

% Create a grid --> change depending on the size of the radius of the
% circle on which the Gaussian optical traps are placed
% r = (4*1e-6) * ones(1,n_particles); % Radial positions for the beams (on circle of radius 12 um)
% theta = linspace(0, (2*pi)-((2*pi)/n_particles), n_particles); % multiple beams evenly spaced around a circle
% [X_translation, Y_translation] = pol2cart(theta, r); % polar coordinate positions of the beams into Cartesian coordinates

% Generate random x and y positions within the FOV
% X_translation = FOV * (rand(1, n_particles) - 0.5);
% Y_translation = FOV * (rand(1, n_particles) - 0.5);

n_particles_per_dim = 3;
% Generate a grid of x and y positions within the FOV
[X_grid, Y_grid] = meshgrid(linspace(-FOV/2+(1/n_particles_per_dim)*FOV, FOV/2-(1/n_particles_per_dim)*FOV, n_particles_per_dim));
% Convert the grid matrices into vectors
X_translation = X_grid(:);
Y_translation = Y_grid(:);

% X_grid = [-50 0 50]*1e-6;
% Y_grid = [-50 0 50]*1e-6;
% [X_grid, Y_grid] = meshgrid(X_grid,Y_grid);
% X_translation = X_grid(:);
% Y_translation = Y_grid(:);

E_beam_individual = zeros(resolution, resolution); % Initialise matrix which will contain the intensities

% Loop through each beam position
for i = 1:numel(X_translation)
    E_beam_individual = E_beam_individual + computeGaussian(lambda_target, w0_target, z_target, [X_translation(i), Y_translation(i)], x_target, y_target);
end

E_target = E_beam_individual;

% Plot the Gaussian beam profile
figure(2);
imagesc(x_target, y_target, abs(E_target).^2);
axis xy;
colormap('default');
title('Target Distribution Beam Profile');
xlabel('X (m)');
ylabel('Y (m)');
colorbar;


%% 3. Implement Gerchberg-Saxton algorithm
% Number of iterations of GS algorithm
n_iterations = 25;

% Initialise E_inv
E_inv = fftshift(ifft2(fftshift(E_target)));

% Initialise error vector
error_vector = zeros(1,n_iterations);

% Convert frequency domain to spatial domain using x = uLf, where f = focal
% length, v = spatial frequency in x direction, L = lambda = wavelength
% f = 2.525e-3; % 2.525 mm
% x_physical = lambda_beam * f * linspace(-0.5,0.5,resolution);
% y_physical = lambda_beam * f * linspace(-0.5,0.5,resolution);

for i = 1:n_iterations

    % Refer to OT textbook page 330 and wikipedia pseudocode for GS
    % algorithm
    E_DOE = abs(E_beam) .* exp(1i*angle(E_inv));
    E_focus = fftshift(fft2(fftshift(E_DOE)));
    E_pre_inv = abs(E_target) .* exp(1i*angle(E_focus));
    E_inv = fftshift(ifft2(fftshift(E_pre_inv)));

    % Normalise the intensities
    max_intensity = max(abs(E_focus(:))); % Find the maximum intensity value
    normalized_E_focus = abs(E_focus) / max_intensity; % Normalize the intensities

    % Present current reconstructed pattern
    figure(3)
    imagesc(x_target,y_target,normalized_E_focus.^2)
    title(['Iteration ' num2str(i)]);
    colorbar
    pause(0.05)

    % Calculate error between target and reconstructed distributions
    error = (sum((abs(E_target(:)) - abs(E_focus(:))).^2))/(resolution*resolution);
    error_vector(i) = error;

end

% Phase mask (phase of CGH) after loop
figure(4);
imagesc(x_slm,y_slm,angle(E_inv))
xlabel('X (m)');
ylabel('Y (m)');
title('Phase mask')
colorbar;

% Reconstruction squared error vs. iterations
iterations = linspace(0,n_iterations,n_iterations);
figure(5)
plot(iterations, error_vector)
title('MSE between target distribution & reconstructed distribution at focal plane')
xlabel('Iterations of GS algorithm')
ylabel('Mean squared error')

%% 4. Combined plots
hold off
tiledlayout(1,2);

% nexttile
% imagesc(x_target, y_target, abs(E_target).^2);
% axis xy;
% colormap('default');
% title('Target Distribution Beam Profile');
% xlabel('X (m)');
% ylabel('Y (m)');
% colorbar;

nexttile
imagesc(x_target,y_target,normalized_E_focus.^2)
title('Final reconstructed focal plane (iteration 25)');
xlabel('X (m)');
ylabel('Y (m)');
colorbar

nexttile
imagesc(x_slm,y_slm,angle(E_inv))
xlabel('X (m)');
ylabel('Y (m)');
title('Phase mask distribution')
colorbar;

%% Function definitions
function intensity = computeGaussian(lambda, w0, z, shift, x, y)

    % Wave number
    k = 2 * pi / lambda;

    % Apply shift (translation of beam in xy plane)
    x = x - shift(1);
    y = y - shift(2);

    % Create 2D mesh with x and y vector
    [X, Y] = meshgrid(x, y);

    % Gaussian beam intensity profile
    wz = w0 * sqrt(1 + (lambda * z / (pi * w0^2))^2);  % Beam waist at distance z
    Rz = z * (1 + (pi * w0^2) / (lambda * z)^2);       % Radius of curvature at distance z
    phi = atan(z / (pi * w0^2) * lambda);              % Gouy phase shift
    
    wz = wz * sqrt(1 + (lambda * z / (pi * wz^2))^2);  % Adjusted beam waist
    intensity = (w0 / wz)^2 * exp(-(X.^2 + Y.^2) / wz^2) * exp(-1i * k * z - 1i * phi);  % Intensity profile

end
