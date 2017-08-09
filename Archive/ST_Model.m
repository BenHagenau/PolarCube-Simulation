%------------------------------
% Model of the PolarCube star tracker
%------------------------------
% Franklin A. Hinckley
% Initial Version: 16 December 2015
% Latest Revision: 1.0
%   16 December 2015
%------------------------------
% Inputs:
%   MRP (3 x 1 double array): true attitude
%   SV (struct): satellite properties
%       .starTrk.noise (double): noise in star tracker output [rad]
% Outputs:
%   MRPsim (3 x 1 double array): simulated star tracker output
%------------------------------

function [MRPsim] = ST_Model(MRP, SV)

% Simulate star tracker output
MRPsim = MRP + SV.starTrk.noise*randn(3,1);
