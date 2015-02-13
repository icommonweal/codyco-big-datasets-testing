function [ output_args ] = hex_calib_to_double( input_args )
    output_args = double(typecast(uint16(sscanf(input_args, '%x')), 'int16'))/double(2^15);
end

