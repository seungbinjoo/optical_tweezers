% md.m - Molecular dynamics simulation of a Brownian particle in a fluid
%
% A microscopic particle (large circle) immersed in a fluid undergoes
% continuous collisions with the fluid molecules (dots). 
%
% See also LJ.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Physical parameters
M = 36;  % number of particles
m = 1;  % particle mass [kg]
e = 1;  % potential depth [J]
rho = .05;  % zero-distance [m]

mb =  100;  % large particle (Brownian particle) mass [kg]
eb = 1;  % potential depth [J]
rhob = .05;  % zero-distance [m]

L = .25;  % box size [m]

%% Simulation parameters
Dt = 5e-5;  % time step [s]
N = 1e+4;  % number of steps

%% Simulation setup

% Position initialization (square lattice)
M = floor(sqrt(M))^2;
r = [ [ceil([.5:1:M]/sqrt(M))]'-.5, [ceil(mod([.5:1:M-.5],sqrt(M)))]'-.5 ]/(sqrt(M))*L; % initial positions [m]

rb = [L/2, L/2];  % Brownian particle initial position [m]

% Momenta inizialization
p = randn(M,2); % [kg m/s]
p = p - ones(M,1)*mean(p);  % remove mean
p = p .* (ones(M,1)*mean(p.^2).^-0.5);  % normalizes variance to 1

pb = [0, 0];

%% Simulation - Leapfrog algorithm
for n = 1:1:N
    % Updates positons (1/2)
    r = r + p/m*Dt/2;

    rb = rb + pb/mb*Dt/2;

    % Claculates forces - Lennard-Jones
    Flj = zeros(M,2);
    for m1 = 1:1:M
        for m2 = m1+1:1:M
            dr = sqrt(sum((r(m2,:)-r(m1,:)).^2));
            Flj(m1,:) = Flj(m1,:) + ( lj(dr+1e-6,e,rho) - lj(dr-1e-6,e,rho) )/2e-6 * [(r(m2,1)-r(m1,1))/dr, (r(m2,2)-r(m1,2))/dr];
            Flj(m2,:) = Flj(m2,:) - ( lj(dr+1e-6,e,rho) - lj(dr-1e-6,e,rho) )/2e-6 * [(r(m2,1)-r(m1,1))/dr, (r(m2,2)-r(m1,2))/dr];
        end
    end

    Fljb = [0, 0];
    for m1 = 1:1:M
        dr = sqrt(sum((rb-r(m1,:)).^2));
            Flj(m1,:) = Flj(m1,:) + ( lj(dr+1e-6,eb,rhob) - lj(dr-1e-6,eb,rhob) )/2e-6 * [(rb(1)-r(m1,1))/dr, (rb(2)-r(m1,2))/dr];
            Fljb = Fljb - ( lj(dr+1e-6,eb,rhob) - lj(dr-1e-6,eb,rhob) )/2e-6 * [(rb(1)-r(m1,1))/dr, (rb(2)-r(m1,2))/dr];
    end

    % Updates momenta
    p = p + Flj*Dt;

    pb = pb + Fljb*Dt;

    % Updates positons (2/2)
    r = r + p/m*Dt/2;

    rb = rb + pb/mb*Dt/2;

    % Reflecting boundaries
    ri = find(r>L | r<0);
    p(ri) = -p(ri);

    rbi = find(rb>L-rhob/2 | rb<rhob/2);
    pb(rbi) = -pb(rbi);

    % Calculates energies
    Vlj = 0;
    for m1 = 1:1:M
        for m2 = m1+1:1:M
            Vlj = Vlj + lj(sqrt(sum((r(m2,:)-r(m1,:)).^2)),e,rho);
        end
        Vlj = Vlj + lj(sqrt(sum((rb-r(m1,:)).^2)),eb,rhob);
    end
    K = sum(sum(0.5*p.^2/m))+sum(0.5*pb.^2/mb);
    E = Vlj+K;

    %% Figure
    if mod(n,10)==0
        
        cla
        title(['Total Energy ' num2str(E) 'J'])
        hold on
        plot(r(:,1),r(:,2),'ko','markerfacecolor','k')
        plot(rb(1),rb(2),'ko','markerfacecolor','r','markersize',40)
        hold on
        axis equal
        axis([-0.05 1.05 -0.05 1.05]*L)
        box on
        set(gca,'XTick',[],'YTick',[])
        drawnow()
        
    end
end

