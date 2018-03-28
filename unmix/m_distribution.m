number_of_fibers = [50 100 250 500 750 1000];
for i = 1:length(number_of_fibers)
    [fibers, fiber_angles] = generate_fibers(number_of_fibers(i));
    cells = generate_cells();
    m = generate_realistic_rt_mixing(fibers, fiber_angles, cells, fiber_profile_exc, fiber_profile_emi, 'figures', false, 'stats', false);
    pd = fitdist(m(:), 'Exponential');
    
    m2 = sort(m, 2, 'descend');
    m2 = m2(:, 1:number_of_fibers(i));
    pd2 = fitdist(m2(:), 'Exponential');
    dist([pd.mu pd2.mu]);
end
