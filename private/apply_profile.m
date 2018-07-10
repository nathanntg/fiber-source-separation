function strength = apply_profile(profile, points)
%APPLY_PROFILE Summary of this function goes here
%   Detailed explanation goes here


if profile.d == 3    
    % calculate strength from sensitivity profile
    strength = interp3(profile.x, profile.y, profile.z, profile.volume, points(1, :), points(2, :), points(3, :));

    % outside
    idx = points(1, :) < profile.x(1) | points(1, :) > profile.x(end) | ...
        points(2, :) < profile.y(1) | points(2, :) > profile.y(end) | ...
        points(3, :) < profile.z(1) | points(3, :) > profile.z(end);

    % clear outside
    strength(idx) = 0;
elseif profile.d == 2
    % convert 3D points to 2D points
    points_r = sqrt(points(1, :) .^ 2 + points(2, :) .^ 2);
    
    % calculate strength from sensitivity profile
    strength = interp2(profile.z, profile.r, profile.volume, points(3, :), points_r);

    % outside
    idx = points_r > profile.r(end) | ...
        points(3, :) < profile.z(1) | points(3, :) > profile.z(end);
    
    % clear outside
    strength(idx) = 0;
else
    error('Unknown profile.');
end

end

