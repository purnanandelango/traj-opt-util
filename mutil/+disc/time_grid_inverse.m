function s = time_grid_inverse(disc,tau,t)
    % Compute the dilation factor, given the time grid and choice of discretization
    % (inverse of time_grid.m)
    N = length(tau);
    switch disc
        case "ZOH"
            N_ctrl = N-1;
        case "FOH"
            N_ctrl = N;
    end
    s = zeros(1,N_ctrl);
    dt = diff(t);
    dtau = diff(tau);
    switch disc
        case "ZOH"
            for k = 1:N_ctrl
                s(k) = dt(k) / dtau(k);
            end
        case "FOH" % NOTE: this is not guaranteed to produce a feasible solution (s \in [s_min,s_max]). For that, solve a separate LP (may be infeasible given the restrictivity of this parameterization).
            s(1) = sum(dt ./ dtau) / length(dt); % boundary condition required to make problem determinant
            for k = 1:N_ctrl-1
                s(k+1) = 2*dt(k)/dtau(k) - s(k);
            end
    end
end