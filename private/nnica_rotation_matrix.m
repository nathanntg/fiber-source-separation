function R = nnica_rotation_matrix(num_sig, axis_i, axis_j, phi)
%NNICA_ROTATION_MATRIX Generate rotation matrix based on axes and angle
%   Generates a rotation matrix that will rotate axis i and j by phi and
%   leave other axes as is. Corresponds with steps 7 in "Algorithms for 
%   Nonnegative Independent Component Analysis" by Mark D. Plumbley.

% leave other axes as is
R = eye(num_sig);

% cosine and sine
phi_cos = cos(phi);
phi_sin = sin(phi);

% update entries
R(axis_i, axis_i) = phi_cos;
R(axis_i, axis_j) = phi_sin;
R(axis_j, axis_i) = - phi_sin;
R(axis_j, axis_j) = phi_cos;

end