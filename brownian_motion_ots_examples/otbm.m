% ot.m - Brownian motion in optical tweezers simulation
%
% Simulates a Brownian motion in a 1D optical trap.
%
% See also BrownianMotion1DOT.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Physical Parameters
dt = 1e-4;
R = 1e-6;
eta = 0.001;
T = 300;
k = 1e-6;
req = 0;

%% Simulation Parameters
N = 1e+5;  % Number of timesteps
r0 = 0;

%% Simulation
bm = BrownianMotion1DOT(dt,R,eta,T,k,req)

bm = bm.simulate(N,r0,'DisplayOn',true)