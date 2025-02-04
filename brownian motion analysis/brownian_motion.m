%% Demonstration of sampling from a Gaussian to produce discrete sequence
%% that simulates white noise --> xyz displacements
%% 1. Simulation of white noise
% Simulation parameters
total_time = 60; % 60 seconds
time_step = 0.1; % time step of 0.1 seconds
num_samples = total_time/time_step;

% Parameters of Gaussian we are sampling from
mean = 0;
variance = 1/time_step;

% Generate Gaussian random numbers with zero mean and unit variance
w_i_x = randn(1,num_samples);
w_i_y = randn(1,num_samples);
w_i_z = randn(1,num_samples);

% Rescale the sequence to obtain variance 1/sqrt(t)
W_i_x = w_i_x / sqrt(time_step);
W_i_y = w_i_y / sqrt(time_step);
W_i_z = w_i_z / sqrt(time_step);

% Plot white noise for xyz directions
figure;
tiledlayout(2,3,"TileSpacing","compact")

nexttile
t = linspace(0,total_time,num_samples);
scatter(t,W_i_x,'MarkerEdgeColor', [0 0.4470 0.7410])
xlabel('t (s)','FontSize',14)
ylabel('W_{n,x}','FontSize',14)

nexttile
scatter(t,W_i_y,'MarkerEdgeColor', [0.9290 0.6940 0.1250])
xlabel('t (s)','FontSize',14)
ylabel('W_{n,y}','FontSize',14)

nexttile
scatter(t,W_i_z,'MarkerEdgeColor', [0.4660 0.6740 0.1880])
xlabel('t (s)','FontSize',14)
ylabel('W_{n,z}','FontSize',14)

%% 2. Resulting particle trajectory in xyz directions due to Brownian motion
r_x = zeros(num_samples,1);
r_y = zeros(num_samples,1);
r_z = zeros(num_samples,1);

for i = 1:num_samples-1
    r_x(i+1) = r_x(i) + sqrt(time_step) * w_i_x(i);
    r_y(i+1) = r_y(i) + sqrt(time_step) * w_i_y(i);
    r_z(i+1) = r_z(i) + sqrt(time_step) * w_i_z(i);
end

% Plot random walk trajectory for xyz directions
nexttile
plot(t,r_x,'Color', [0 0.4470 0.7410])
xlabel('t (s)','FontSize',14)
ylabel('x_{n}','FontSize',14)

nexttile
plot(t,r_y,'Color', [0.9290 0.6940 0.1250])
xlabel('t (s)','FontSize',14)
ylabel('y_{n}','FontSize',14)

nexttile
plot(t,r_z,'Color', [0.4660 0.6740 0.1880])
xlabel('t (s)','FontSize',14)
ylabel('z_{n}','FontSize',14)