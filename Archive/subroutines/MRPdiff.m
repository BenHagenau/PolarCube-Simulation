%------------------------------
% Computes the relative orientation of two frames in MRPs
% Difference is frame 2 relative to frame 1
%------------------------------
% Franklin Hinckley
% 28 August 2015
%------------------------------
% Inputs:
%   MRP1 (3 x 1 double matrix): Modified Rodrigues Parameters
%   MRP2 (3 x 1 double matrix): Modified Rodrigues Parameters
% Outputs:
%   MRP (3 x 1 double matrix): difference MRP
%------------------------------

function [MRP] = MRPdiff(MRP1, MRP2)

% Compute DCM for MRP1
DCM1 = mrp2dcm(MRP1);

% Compute DCM for MRP2
DCM2 = mrp2dcm(MRP2);

% Compute difference between frames
DCM = DCM2 * DCM1';

% Extract quaternion
[ep] = dcm2ep(DCM);

% Convert to MRP
MRP = ep2mrp(ep);

end