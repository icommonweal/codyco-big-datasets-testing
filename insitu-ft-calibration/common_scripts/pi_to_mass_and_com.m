function [ mass_and_com ] = pi_to_mass_and_com( pi )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    mass_and_com = pi;
    mass_and_com(2:4) = pi(2:4)./pi(1);

end

