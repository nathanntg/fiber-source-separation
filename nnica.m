function [Z, M, W] = nnica(Y)
%NNICA Summary of this function goes here
%   Detailed explanation goes here

% number of dimensions to extract
num_sig = size(Y, 1);
max_iter = 100000;
tolerance = 1e-5;
verbose = 'on';

b_verbose = strcmp(verbose, 'on');

% calculate PCA
[E, D] = pcamat(Y, 1, num_sig, 'off', verbose);

% Calculate the whitening
[Y_white, whiten, dewhiten] = whitenv(Y, E, D, verbose);

% in case number of signals was decreased during whitening
num_sig = size(Y_white, 1);

% initial W matrix
W = eye(num_sig); % unmixing matrix
Z = Y_white; % unmixed data

% torque matrix
T = nnica_torque_matrix(Z);

% loop
torque = nan;
for i = 1:max_iter
    %Z = W * Z;
    
    % figure out maximal torque, and accompanying row / column
    [torque, idx] = max(T(:));
    [axis_i, axis_j] = ind2sub(size(T), idx);
    if torque < tolerance
        break;
    end
    
    % find angle
    phi = nnica_min_rotation_error(Z([axis_i, axis_j], :));
    
    % perform rotation
    R = nnica_rotation_matrix(num_sig, axis_i, axis_j, phi);
    
    % rotate data
    W = R * W;
    Z = R * Z;
    
    % update torques
    T = nnica_torque_update(Z, T, [axis_i, axis_j]);
end

% combining whitening and unmixing matrices
% and generate mixing matrix
M = dewhiten * W';
W = W * whiten;

% print stopping conditions
if torque < tolerance
    if b_verbose
        fprintf('Stopped at iteration %d with torque %f, below tolerance.\n', i, torque);
    end
else
    warning('Solution did not converge after %d iterations.', max_iter);
end

end
