function [R_hat, theta_hat, v_hat] = main_MIMO(num_PRIs, SNR_dB, target_range, target_angles_deg, target_velocity, virtual, Ntx, Nrx)

    A = 150;
    Tc = 28e-9;
    M_range_bins = 199;
    Nc = 127;
    c = 3e8;
    Fc = 15e9;
    lambda = c / Fc;
    d = lambda / 2;

    % variable array geometry
    x_rx = ((0:Nrx-1).' - (Nrx-1)/2) * d;
    x_tx = ((0:Ntx-1).' - (Ntx-1)/2) * Nrx * d;
    x = x_rx;
    gain = Ntx;

    target_rcs = [1 1];
    num_targets = length(target_range);
    angle_range = 30:0.1:150;

    [baseband_signal, mseq_codes] = gen_Tx_signal_MIMO(A, num_PRIs, M_range_bins, Ntx);

    backscatter_data = gen_backscatter_MIMO(baseband_signal, lambda, ...
                                    M_range_bins, Nc, num_PRIs, target_range, ...
                                    target_angles_deg, x, x_tx, x_rx, Tc, ...
                                    target_rcs, target_velocity, gain, Fc, virtual);

    [noise, noise_power] = gen_noise_MIMO(backscatter_data, SNR_dB);

    backscatter_data = backscatter_data + noise;

    [R_hat, theta_hat, v_hat] = find_target_MIMO( ...
                                    angle_range, backscatter_data, mseq_codes, ...
                                    M_range_bins, Nc, num_PRIs, ...
                                    x, x_tx, x_rx, lambda, Tc, num_targets, virtual, ...
                                    noise_power);
end