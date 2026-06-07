% Generate backscatter data from Tx signal
function backscatter_data = gen_backscatter_MIMO(Tx_signal, lambda, ...
                                    M_range_bins, Nc, num_PRIs, target_range, ...
                                    target_angles_deg, x, x_tx, x_rx, Tc, ...
                                    target_rcs, target_velocity, gain, Fc, virtual)

    c = 3e8;
    PRI_len = M_range_bins * Nc;
    T_PRI = PRI_len * Tc;
    total_samples = PRI_len * num_PRIs;
    Nrx = length(x_rx);
    Ntx = length(x_tx);
    backscatter_data = zeros(Nrx, total_samples);

    for i = 1:length(target_angles_deg)

        if virtual
            phi_tx = ((2*pi) / lambda) * x_tx(:) * cosd(target_angles_deg(i));
            phi_rx = ((2*pi) / lambda) * x_rx(:) * cosd(target_angles_deg(i));
    
            manifold_tx = exp(+1j * phi_tx);
            manifold_rx = exp(-1j * phi_rx);
        else
            phi_target = ((2*pi) / lambda) * x * cosd(target_angles_deg(i));
            manifold_tx = exp(+1j * phi_target);
            manifold_rx = exp(-1j * phi_target);
        end

        % transmit projection towards target
        X_tx = manifold_tx' * Tx_signal; 

        % target delay
        t_echo = (2 * target_range(i)) / c;
        delay_samples = round(t_echo / Tc);

        delayed_signal = zeros(1, total_samples);

        if delay_samples < total_samples
            delayed_signal(1 + delay_samples : total_samples) = ...
                X_tx(1 : total_samples - delay_samples);
        end

        % doppler
        delayed_pri = reshape(delayed_signal, PRI_len, num_PRIs);

        fd = 2 * target_velocity(i) / lambda;

        for p = 1:num_PRIs
            doppler_phase = exp(1j * 2*pi * fd * (p-1) * T_PRI);
            delayed_pri(:,p) = delayed_pri(:,p) * doppler_phase;
        end

        delayed_signal = reshape(delayed_pri, 1, []);

        % RCS model
        if i == 1
            rcs = target_rcs(i);
        elseif i == 2
            sigma_mean = target_rcs(i);
            B = sqrt(sigma_mean / 2);
            amp = raylrnd(B, [1 total_samples]);
            rcs = amp.^2;
        elseif i == 3
            rcs = target_rcs(i)/4 * chi2rnd(4, [1 total_samples]);
        end

        % scattering coefficient
        random_phase = 2*pi*rand(1,1);

        beta_1 = sqrt(gain.^2 / (4*pi)^3) ...
                 * (lambda / target_range(i)^2) ...
                 * sqrt(rcs);

        beta = beta_1 ...
               * exp(-1j * 2*pi * Fc * t_echo) ...
               .* exp(1j * random_phase);

        signal = beta .* delayed_signal;

        % receive array manifold
        signal_w_mani = manifold_rx * signal;   
        backscatter_data = backscatter_data + signal_w_mani;
    end
end