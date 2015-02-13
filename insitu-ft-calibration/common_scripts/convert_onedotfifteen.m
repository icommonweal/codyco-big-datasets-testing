%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
% Function to convert a decimal number into 1.15 hex format
function [hex_out] = convert_onedotfifteen(indecimal)

temp = indecimal; % Ensure positive number
quant = abs(round((2^15)*-temp)); % Quantize to 15 bits and round to the nearest integer;

if indecimal < 0 % i.e. negative
quant = bitcmp(quant,16) + 1; % Take the 2's complement
end

hex_out=dec2hex(quant); % Convert the decimal number to hex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%