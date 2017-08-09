%--------------------------------
% Computes the product of two MRPs
%   Result is MRP1 x MRP2
%--------------------------------
% Franklin Hinckley
% 16 December 2015
%--------------------------------
% Inputs:
%   MRP1 (3 x 1 double array): first MRP
%   MRP2 (3 x 1 double array): second MRP
% Outputs:
%   MRP (3 x 1 double array): MRP product (composite rotation)
%--------------------------------

function [MRP] = MRPmult(MRP1, MRP2)

% Multiply MRPs
MRP = ((1-MRP2'*MRP2)*MRP1 + (1-MRP1'*MRP1)*MRP2 - 2*tilde(MRP1)*MRP2)/...
    (1 + (MRP2'*MRP2)*(MRP1'*MRP1) - 2*MRP1'*MRP2);

end
