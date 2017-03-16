function profile = sp_3d_to_2d(profile)
%SP_3D_TO_2D Convert 3D profile to 2D profile
%   If a 3D profile is symmetric around the Z asxis, then it can be
%   converted to a 2D profile for faster processing.

% check input
if profile.d ~= 3
    error('Expected a valid 3D profile.');
end

% find center
idx_zero = find(profile.y == 0); % zero Y value
idx_positive = profile.x >= 0; % positive X values

% ensure a single center point was found
if ~isscalar(idx_zero)
    error('Could not find y = 0 axis for slicing profile.');
end

% create new profile
new_volume = squeeze(profile.volume(idx_zero, idx_positive, :));

% create new profile
profile = struct('r', profile.x(idx_positive), ...
    'z', profile.z, 'volume', new_volume, 'd', 2);

end
