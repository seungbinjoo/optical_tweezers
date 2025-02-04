% Parameters
lambda = 1e-6;          % Wavelength (in meters)
k = 2*pi/lambda;        % Wave number
w0 = 1e-3;              % Waist radius of Gaussian beam (in meters)
z = 0;                  % Propagation distance (in meters)
Nx = 100;               % Number of points in x direction
Ny = 100;               % Number of points in y direction
Lx = 5*w0;              % Range of x-axis (in meters)
Ly = 5*w0;              % Range of y-axis (in meters)

% Generate meshgrid
x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(-Ly/2, Ly/2, Ny);
[X, Y] = meshgrid(x, y);

% Gaussian Beam
wz = w0 * sqrt(1 + (lambda*z/(pi*w0^2))^2);   % Beam radius at distance z
Rz = z * (1 + (pi*w0^2)/(lambda*z)^2);         % Radius of curvature at distance z
phase = k*(X.^2 + Y.^2)/(2*Rz) - atan(z/lambda); % Gouy phase
Gaussian_beam = exp(-(X.^2 + Y.^2)/wz^2) .* exp(-1j * phase);

% Laguerre-Gaussian Beam (p = 1, l = 0)
p = 1; l = 0;
Laguerre_Gaussian_beam = sqrt(2/pi) * sqrt(factorial(l)/(pi*factorial(l+p))) .* ...
    (sqrt(2) * sqrt(X.^2 + Y.^2)/w0).^l .* exp(-((X.^2 + Y.^2)/(w0^2))) .* ...
    Lpoly(p, 2*l, 2*((X.^2 + Y.^2)/w0^2));

% Hermite-Gaussian Beam (n = 0, m = 1)
n = 0; m = 1;
Hermite_Gaussian_beam = hermite_poly(n, sqrt(2)*X/w0) .* hermite_poly(m, sqrt(2)*Y/w0) ...
    .* exp(-(X.^2 + Y.^2)/w0^2);

% Plotting
figure;
subplot(1,3,1);
imagesc(x, y, abs(Gaussian_beam).^2);
title('Gaussian Beam');
xlabel('x (m)');
ylabel('y (m)');
axis square;

subplot(1,3,2);
imagesc(x, y, abs(Laguerre_Gaussian_beam).^2);
title('Laguerre-Gaussian Beam (p=1, l=0)');
xlabel('x (m)');
ylabel('y (m)');
axis square;

subplot(1,3,3);
imagesc(x, y, abs(Hermite_Gaussian_beam).^2);
title('Hermite-Gaussian Beam (n=0, m=1)');
xlabel('x (m)');
ylabel('y (m)');
axis square;

function L = Lpoly(p, l, r)
    L = zeros(size(r));
    for k = 0:l
        L = L + (-1)^k * nchoosek(l,k) * factorial(l+p) / (factorial(k) * factorial(p+k) * factorial(l-k) * factorial(p)) * r.^k;
    end
    L = L .* (sqrt(2) * r).^p .* exp(-r.^2/2);
end

function H = hermite_poly(n, x)
    H = (-1)^n * exp(x.^2) .* diff(exp(-x.^2), n);
end