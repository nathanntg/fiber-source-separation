function m = generate_mixing_matrix(inputs, outputs)
%GENERATE_MIXING_MATRIX Exponentially distributed random mixing matrix
%   The most naive model, just treats the mixing matrix as the result of an
%   exponential distribution.

m = exprnd(1.5, [outputs, inputs]);

% normalize (can't "increase" the amount of light)
m = m / max(max(m));

end

