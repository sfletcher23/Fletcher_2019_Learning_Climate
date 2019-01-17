function [out] = cmps2cmpd(in)

% in is a matrix or scalar in cubic meters per sec
% out is a matrix of same size in cubmic meters per day

out = in * 60 *60 * 24;
