function c = generate_smoothing(rows, columns, gen)
%GENERATE_SMOOTHING Summary of this function goes here
%   Detailed explanation goes here

l_rows = length(rows);
l_columns = length(columns);
if l_rows == 1
    num = l_columns;
    rows = ones(1, num) * rows;
elseif l_columns == 1
    num = l_rows;
    columns = ones(1, num) * columns;
elseif l_rows == l_columns
    num = l_rows;
else
    error('Dimensions for smoothing matrix must be the same, or either must be 1');
end

if nargin < 3
    gen = @(n, m) fspecial('average', [n m]); %normrnd(0.5, 0.25, n, m);
end

% make cell array
c = cell(1, num);
for i = 1:num
    c{i} = gen(rows(i), columns(i));
end

end

