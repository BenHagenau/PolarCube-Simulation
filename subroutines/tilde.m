%------------------------------
% Computes the tilde matrix for the provided vector
%------------------------------
% Franklin Hinckley
% 30 October 2015
%------------------------------
% Inputs:
%   vec (3x1 double vector): input vector
% Outputs:
%   tilde (3x3 double matrix): tilde matrix for the vector
%------------------------------

function [tilde]=tilde(vec)

%Compute the tilde matrix
tilde=[0 -vec(3) vec(2);...
    vec(3) 0 -vec(1);...
    -vec(2) vec(1) 0];