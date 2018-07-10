function fibers_max = max_profile_intensity(fibers, fiber_angles, profile)
%MAX_PROFILE_INTENSITY Summary of this function goes here
%   Detailed explanation goes here

[~, ind] = max(profile.volume(:));
if profile.d == 2
    [r, z] = ind2sub(size(profile.volume), ind);
    rel = [0; profile.r(r); profile.z(z)];
else
    [x, y, z] = ind2sub(size(profile.volume), ind);
    rel = [profile.x(x); profile.y(y); profile.z(z)];
end

fibers_max = zeros(size(fibers));
for i = 1:size(fibers, 2)
    theta = fiber_angles(:, i); % angle
    r = [1 0 0; 0 cos(theta(1)) -sin(theta(1)); 0 sin(theta(1)) cos(theta(1))]; % x
    r = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))] * r; % y
    
    fibers_max(:, i) = fibers(:, i) + r * rel;
end

end

