% stimulation threshold
light_out = 0.0544; % uW out of the fiber
threshold = 5; % mW / mm^2

% above stimulation threshold
stim = profile_exc.volume * light_out  > (threshold / 1000);

% convert to cubic microns
volume = sum(stim(:));
volume = volume * median(diff(profile_exc.x)) * median(diff(profile_exc.y)) * median(diff(profile_exc.z));
fprintf('%d um^3\n', volume);

% make slice
idx_zero = find(profile_exc.y == 0); % zero Y value
slice = squeeze(stim(idx_zero, :, :));

imagesc(profile_exc.z, profile_exc.x, slice)
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');
title(sprintf('Stimulation threshold (%d mW / mm^2)', threshold));
