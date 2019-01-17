function [out] = cmpd2mcmpy(in)

% in is a matrix or scalar in cubic meters per day
% out is a matrix of same size in million cubmic meters per year

out = in * 365 / 1E6;