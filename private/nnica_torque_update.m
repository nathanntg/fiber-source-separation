function T = nnica_torque_update(Z, T, axes_to_update)
%NNICA_TORQUE_UPDATE Update a torque matrix
%   Updates an existing torque matrix, used to recalculate the torque
%   after a rotation is applied to a set of axes. Expects an existing
%   torque matrix generated by NNICA_TORQUE_MATRIX, where only half of the
%   matrix is populated (where row < col). Torques are recalculated for all
%   axes that include one of axes_to_update. The torque matrix can be used 
%   to identify dimension pairs that produce maximal torque. Corresponds
%   with steps 2-3 in "Algorithms for Nonnegative Independent Component 
%   Analysis" by Mark D. Plumbley.

% number of signals
num_sig = size(Z, 1);

% rectified portions
Zp = max(Z, 0);
Zn = min(Z, 0);

% iterate over axes needing update
for k = 1:length(axes_to_update)
    i = axes_to_update(k);
    
    % iterate over pairs
    for j = 1:num_sig
        % no calculation needed if same index
        if i == j
            continue;
        end
        
        % calculate torque
        torque = abs(sum(Zp(i, :) .* Zn(j, :) - Zn(i, :) .* Zp(j, :)));
        
        % keep consistent half of the torque matrix filled
        if i < j
            T(i, j) = torque;
        else
            T(j, i) = torque;
        end
    end
end

end
