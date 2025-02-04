%% Workspace initialisation

clear all;
close all;
clc;

%% Plot autocorrelation
k_B = 1.380e-23; % Boltzmann constant
T = 293; % temperature (K)
kappa = [0.1 0.5 1 5 10] * (1e-15/1e-9); % trap stiffnesses (fN/nm)

R = 1e-6; % Particle radius [m] --> 1um radius particle
eta = 0.001; % Water viscosity [Pa*s]
gamma = (6*pi*eta*R);
tau_to = gamma./kappa;

for i = 1:numel(kappa)
    tau = -500e-3:1e-3:500e-3;
    acf = ((k_B*T)/(kappa(i)))*exp(-abs(tau)/tau_to(i));
    acf_normalized = acf / max(acf);  % Normalize to maximum value of 1
    hold on
    plot(tau,acf_normalized,'LineWidth',1)
    xlim([-500e-3 500e-3])
    ylim([0 1])
    grid
    xlabel('\tau (s)')
    ylabel('Autocorrelation R_{xx}')
    title('Position autocorrelation function of optically trapped Brownian particle')
end

legend('1 fN/nm', '5 fN/nm', '10 fN/nm', '15 fN/nm', '20 fN/nm')