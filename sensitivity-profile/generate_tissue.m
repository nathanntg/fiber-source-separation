function [ output_args ] = generate_tissue(fn, nx, ny, nz)
%GENERATE_TISSUE Summary of this function goes here
%   Detailed explanation goes here

% load second 2pt
fid = fopen(fn, 'wb');
fwrite(fid, zeros(1, nx * ny, 'uint8'), 'uint8');
fwrite(fid, ones(1, nx * ny * (nz - 1), 'uint8'), 'uint8');
fclose(fid);


end

