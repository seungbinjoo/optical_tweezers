%% Workspace initialisation
clear all;
close all;
clc;

%% Plot MSD
k_B = 1.380e-23; % Boltzmann constant
T = 293; % temperature (K)
kappa = [0.2 0.83 3 5.83 15] * (1e-15/1e-9); % trap stiffnesses (fN/nm)

R = 1e-6; % Particle radius [m] --> 1um radius particle
eta = 0.001; % Water viscosity [Pa*s]
gamma = (6*pi*eta*R);
tau_to = gamma./kappa;
msd = zeros(numel(kappa),50);

figure; % Create a new figure window
for i = 1:numel(kappa)
    tau = logspace(-3,0);
    msd(i,:) = (2*k_B*T)/(kappa(i))*(1-exp(-abs(tau)./tau_to(i)));
end

msd = msd * 10^(12); % conversion from m^2 to um^2

loglog(tau,msd(1,:),tau,msd(2,:),tau,msd(3,:),tau,msd(4,:),tau,msd(5,:),'LineWidth',1)
grid on
xlabel('\tau (s)')
ylabel('MSD (\mum^2)')
title('Mean square displacement of optically trapped Brownian particle')

legend('0.2 pN/\mum', '0.83 pN/\mum', '3 pN/\mum', '5.83 pN/\mum', '15 pN/\mum','Location', 'northwest')