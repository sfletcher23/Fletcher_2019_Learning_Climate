function [out] = mcmpy2cmpd(in)

% in is a matrix or scalar in million cubmic meters per year
% out is a matrix of same size in cubic meters per day

out = in * 1E6 / 365;