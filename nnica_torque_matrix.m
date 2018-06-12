function T = nnica_torque_matrix(Z)
%NNICA_TORQUE_MATRIX Build a torque matrix
%   Calculates the torque between dimension pairs in matrix Z. The
%   resulting torque matrix is symmetric, but this function only populates 
%   one half of the matrix (in theory, this data could be stored in a 
%   vector form like pdist, but probably not worth while to keep moving 
%   between vector and matrix coordinates). The torque matrix can be used 
%   to identify dimension pairs that produce maximal torque. Corresponds
%   with steps 2-3 in "Algorithms for Nonnegative Independent Component 
%   Analysis" by Mark D. Plumbley.

% number of signals
num_sig = size(Z, 1);

% rectified portions
Zp = max(Z, 0);
Zn = min(Z, 0);

% empty torque matrix
T = zeros(num_sig, num_sig);

% iterate over axis pairs
for i = 1:(num_sig - 1)
    for j = (i + 1):num_sig
        T(i, j) = abs(sum(Zp(i, :) .* Zn(j, :) - Zn(i, :) .* Zp(j, :)));
    end
end

end
