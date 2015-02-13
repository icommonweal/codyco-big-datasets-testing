function [ output_args ] = denormalize2( x,y,z )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
     output_args = (x.*(ones(size(x,1),1)*z)+ones(size(x,1),1)*y);


end

