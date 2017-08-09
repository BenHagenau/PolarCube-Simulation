%--------------------------------
% Use Symbolic Toolbox to determine STM and PNTM
%--------------------------------
% Franklin Hinckley
% 23 December 2015
%--------------------------------

%% Setup
clearvars
close all
clc

addpath subroutines

%% Define variables
s1 = sym('s1','real');
s2 = sym('s2','real');
s3 = sym('s3','real');

ws1 = sym('ws1','real');
ws2 = sym('ws2','real');
ws3 = sym('ws3','real');

b1 = sym('b1','real');
b2 = sym('b2','real');
b3 = sym('b3','real');

dt = sym('dt','real');

sigma = [s1 s2 s3]';
wStar = [ws1 ws2 ws3]';
b = [b1 b2 b3]';

%% Define kinematics
%B = ((1-(sigma'*sigma)))*eye(3)+(2*tilde(sigma))+(2*(sigma*sigma'));

%% Define f,g vectors
f = [0.25*(((1-(sigma'*sigma)))*eye(3)+(2*tilde(sigma))+(2*(sigma*sigma')))*(wStar-b); zeros(3,1)];
%g = [-0.25*B(sigma)*etaW; etaB];

%% Compute partials
A = jacobian(f,[sigma; b]);

%% Compute STM
Phi = exp(A*dt);
