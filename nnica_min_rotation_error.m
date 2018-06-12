function phi = nnica_min_rotation_error(Zs)
%NNICA_MIN_ROTATION_ERROR Line search on rotation to minimize error
%   Find the angle phi that minimizes the error term associated with 
%   rotating the two dimensional data into the positive quadrant. 
%   Corresponds with steps 6 in "Algorithms for Nonnegative Independent
%   Component Analysis" by Mark D. Plumbley.

    % check size
    assert(size(Zs, 1) == 2);
    
%     % line search
%     phi1 = line_search(64);
%     phi2 = quadrant_search();
%     phi3 = matlab_search();
%     
%     disp([phi1 phi2 phi3]);
%     disp([test(phi1) test(phi2) test(phi3)]);
%     
%     phi = phi3;
    
    phi = matlab_search();
    
    function phi = line_search(steps)
        phi = 0;
        best = inf;
        for c_phi = linspace(0, 2 * pi, steps)
            % rotation matrix
            c_phi_cos = cos(c_phi);
            c_phi_sin = sin(c_phi);
            W = [c_phi_cos c_phi_sin; - c_phi_sin c_phi_cos];
            
            % apply rotation
            Y = W * Zs;
            
            % calculate error
            err = J(W, Y);
            
            % update best error
            if err < best
                phi = c_phi;
                best = err;
            end
        end
    end

    function c_phi = quadrant_search()
        tol = 1e-4;
        maxiter = 10;
        c_phi = 0;
        for i = 1:maxiter
            % rotation matrix
            c_phi_cos = cos(c_phi);
            c_phi_sin = sin(c_phi);
            W = [c_phi_cos c_phi_sin; - c_phi_sin c_phi_cos];
            
            % apply rotation
            Y = W * Zs;
            
            % update
            c_phi_new = c_phi - 2 * J(W, Y) / dJ(Y);
            delta = abs(c_phi_new - c_phi);
            c_phi = c_phi_new;
            
            if delta < tol
                break;
            end
        end
    end

    function phi = matlab_search()
        phi = fzero(@phi2dJ, line_search(8));
    end

    function dJv = phi2dJ(c_phi)
        % rotation matrix
        c_phi_cos = cos(c_phi);
        c_phi_sin = sin(c_phi);
        W = [c_phi_cos c_phi_sin; - c_phi_sin c_phi_cos];

        % apply rotation
        Y = W * Zs;
        
        dJv = dJ(Y);
    end

    function v = test(c_phi)
        % rotation matrix
        c_phi_cos = cos(c_phi);
        c_phi_sin = sin(c_phi);
        W = [c_phi_cos c_phi_sin; - c_phi_sin c_phi_cos];
        Y = W * Zs;
        
        v = [J(W, Y); dJ(Y)];
    end

    function v = J(W, Y)
        % square error
        e2 = (Zs - W' * max(Y, 0)) .^ 2;
        v = 0.5 * mean(e2(:));
    end

    function v = dJ(Y)
        % rectify data
        Yp = max(Y, 0);
        Yn = min(Y, 0);
        
        v = -mean(Yp(1, :) .* Yn(2, :) - Yn(1, :) .* Yp(2, :));
    end
end

