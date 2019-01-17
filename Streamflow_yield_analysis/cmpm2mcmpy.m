function [out] = cmpm2mcmpy(in)

% in is a matrix or scalar in cubic meters per month
% out is a matrix of same size in million cubmic meters per year

out = in * 12 / 1E6;