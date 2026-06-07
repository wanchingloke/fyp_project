function [R_hat, theta_hat, v_hat] = find_target_MIMO( ...
                                            angle_range, backscatter_data, mseq_codes, ...
                                            M_range_bins, Nc, num_PRIs, ...
                                            x, x_tx, x_rx, lambda, Tc, num_targets, virtual_mode, ...
                                            noise_power)

    c = 3e8;
    range_res = c * Tc / 2;
    [Nrx, total_samples] = size(backscatter_data);
    [Ntx, code_len] = size(mseq_codes);
    Nc = code_len;
    PRI_len = M_range_bins * Nc;

    % reshape received signal into radar data cube
    rx_cube = reshape(backscatter_data, Nrx, PRI_len, num_PRIs);

    % matched filtering with each PN code
    Y = zeros(Nrx, PRI_len, num_PRIs, Ntx);

    for tx = 1:Ntx
        code = mseq_codes(tx,:);
        h = conj(fliplr(code));
        for p = 1:num_PRIs
            for rx = 1:Nrx
                Y(rx,:,p,tx) = conv(rx_cube(rx,:,p), h, 'same') / Nc;
            end
        end
    end

    % detection thresholding
    P_range = squeeze(sum(sum(sum(abs(Y).^2, 1), 3), 4));
    P_range = P_range(:);
    decision_threshold = (raylinv(0.99, sqrt(noise_power)))^2;
    detection = P_range > decision_threshold;
    P_cfar = P_range .* detection;
  
    % find detected range bins
    min_range_sep_m = 20;
    min_range_sep_bins = round(min_range_sep_m / range_res);
    [~, pkBins] = findpeaks(P_cfar, ...
                        'SortStr','descend', ...
                        'NPeaks',num_targets, ...
                        'MinPeakDistance',min_range_sep_bins);
    
    if isempty(pkBins)
        fprintf("No range peaks found using CFAR.\n");
        R_hat = [];
        theta_hat = [];
        v_hat = [];
        return;
    end
    
    % initialise outputs
    R_hat = [];
    theta_hat = [];
    v_hat = [];
   
    T_PRI = PRI_len * Tc;
    num_range_bins_detected = length(pkBins);
    
    % process each detected range bin
    for k = 1:length(pkBins)
    
        rbin = pkBins(k);
    
        % range estimate
        delay_samples = rbin - floor(Nc/2) - 1;
        R_current = delay_samples * range_res;
    
        % doppler estimate
        slow_signal = squeeze(sum(sum(Y(:,rbin,:,:), 1), 4));
        slow_signal = slow_signal(:);
        Sdop = fftshift(fft(slow_signal, num_PRIs));
        fd_axis = (-num_PRIs/2 : num_PRIs/2-1).' / (num_PRIs * T_PRI);
        [~, dop_idx] = max(abs(Sdop));
        fd_hat = fd_axis(dop_idx);
        v_current = lambda * fd_hat / 2;
        
        phase_comp = exp(-1j * 2*pi * fd_hat * (0:num_PRIs-1) * T_PRI);

        if virtual_mode == false
            % classical MIMO mode

            Z = zeros(Nrx, Ntx);
            phase_comp = exp(-1j * 2*pi * fd_hat * (0:num_PRIs-1) * T_PRI);

            for tx = 1:Ntx
                snaps = squeeze(Y(:,rbin,:,tx));   % Nrx x num_PRIs
                snaps_comp = snaps .* phase_comp;
                Z(:,tx) = mean(snaps_comp, 2);
            end

            % classical beamforming
            score = zeros(size(angle_range));

            for a = 1:length(angle_range)

                theta = angle_range(a);
                phi = (2*pi/lambda) * x(:) * cosd(theta);
                S_R = exp(-1j * phi);
                S_T = exp(+1j * phi);
                beta_LS(a) = (S_R' * Z * S_T) / (Nrx * Ntx);
                score(a) = abs(beta_LS(a));

            end
        else
            % virtual MIMO mode

            Z = zeros(Nrx, Ntx);
            for tx = 1:Ntx
                snaps = squeeze(Y(:,rbin,:,tx));   
                snaps_comp = snaps .* phase_comp;
                Z(:,tx) = mean(snaps_comp, 2);
            end

            % virtual beamforming
            z_v = Z(:);
            score = zeros(size(angle_range));
            
            for a = 1:length(angle_range)
            
                theta = angle_range(a);
                phi_rx = (2*pi/lambda) * x_rx(:) * cosd(theta);
                phi_tx = (2*pi/lambda) * x_tx(:) * cosd(theta);
                S_R = exp(-1j * phi_rx);
                S_T = exp(+1j * phi_tx);
                S_v = kron(conj(S_T), S_R);
                score(a) = abs(S_v' * z_v)^2 / (Nrx*Ntx)^2;
            
            end
        end

        % detect angle peaks
        angle_threshold = 0.15 * max(score);
        
        [angle_pks, angle_locs] = findpeaks(score, ...
            'SortStr','descend', ...
            'MinPeakDistance',2, ...
            'MinPeakHeight',angle_threshold);
        
        if isempty(angle_locs)
            [~, angle_locs] = max(score);
        end
        
        % limit detections
        remaining_targets = num_targets - length(R_hat);
        
        if length(angle_locs) > remaining_targets
            angle_locs = angle_locs(1:remaining_targets);
        end
        
        % store detections
        for aa = 1:length(angle_locs)
        
            theta_current = angle_range(angle_locs(aa));
        
            R_hat(end+1,1) = R_current;
            theta_hat(end+1,1) = theta_current;
            v_hat(end+1,1) = v_current;
        
            fprintf("Target %d: Range = %.2f m , Angle = %.2f deg , Velocity = %.2f m/s\n", ...
                length(R_hat), R_current, theta_current, v_current);

        end

    end
end