function m = generate_mixing_matrix(inputs, outputs)
%GENERATE_MIXING_MATRIX Summary of this function goes here
%   Detailed explanation goes here

m = exprnd(1.5, [outputs, inputs]);

% normalize (can't "increase" the amount of light)
m = m / max(max(m));

end

