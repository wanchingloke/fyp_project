function r_hat = laa_single(SNR_dB, target, r, N, alpha)

    % constants
    Fc = 2.4e9;         % carrier freq = 2.4GHz
    c = 3e8;            % velocity of light
    lambda = c / Fc;    % wavelength
    k = 2*pi/lambda;    % wavenumber
    N_rx = size(r,2);

    % single target signal
    n = 0:N-1;
    s = exp(1j*2*pi*0.05*n);

    % true ranges
    rho = zeros(N_rx,1);
    for i = 1:N_rx
        rho(i) = norm(target - r(:,i));
    end

    gamma = zeros(N_rx,1);
    noise_var = zeros(N_rx,1);

    % rotate LAA reference point
    for ref = 1:N_rx

        rho_ref = rho(ref);

        % spherical manifold vector
        a = (rho_ref ./ rho).^alpha .* ...
            exp(-1j*k*(rho - rho_ref));

        % clean received signal
        x_clean = a * s;

        % add noise
        signal_power = mean(abs(x_clean(:)).^2);
        noise_power = signal_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power/2) * (randn(size(x_clean)) + 1j*randn(size(x_clean)));
        x = x_clean + noise;

        % covariance matrix
        Rxx = (x * x') / N;

        % eigenvalue decomposition
        eigvals = sort(real(eig(Rxx)), 'descend');
        gamma(ref) = eigvals(1);

        % everything but first eigenvector is noise subspace
        noise_var(ref) = mean(eigvals(2:end));

    end

    % LAA signal eigenvalues
    lambda_sig = gamma - noise_var;
    lambda_sig(lambda_sig <= 0) = eps;

    % range-ratio K
    K_hat = (lambda_sig(2:end) / lambda_sig(1)).^(1/(2*alpha));

    % metric fusion stage
    r_hat = laa_metric_fusion(r, K_hat);

    % print output
    fprintf('True position: x = %.2f m, y = %.2f m\n', target(1), target(2));
    fprintf('Estimated position: x = %.2f m, y = %.2f m\n', r_hat(1), r_hat(2));
    fprintf('Position error: %.3f m\n', norm(r_hat - target));

end


function r_hat = laa_metric_fusion(r, K)

    N_rx = size(r,2);
    r_ref = r(:,1);
    r_rest = r(:,2:end);

    H = zeros(N_rx-1, 3);
    b = zeros(N_rx-1, 1);

    for i = 1:N_rx-1

        ri = r_rest(:,i);
        Ki = K(i);

        H(i,:) = [2*(r_ref - ri).', 1 - Ki^2];

        b(i) = norm(r_ref)^2 - norm(ri)^2;

    end
    
    % solve metric fusion
    z_hat = H \ b;
    r_hat = z_hat(1:2);

end