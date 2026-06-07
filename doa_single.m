function [r_hat, theta_hat, rho_hat] = doa_single(SNR_dB, target, n_elements, N_snapshots, r)

    % constants
    Fc = 2.4e9;         % carrier freq = 2.4GHz
    c = 3e8;            % velocity of light
    lambda = c / Fc;    % wavelength
    k = 2*pi/lambda;    % wavenumber
    N_rx = size(r,2);   % number of receivers

    % UCA antenna positions (n_elements antennas are placed evenly around a
    % circle)
    radius = 0.125; 
    theta = linspace(0, 2*pi, n_elements+1);
    theta(end) = [];
    rUCA = [radius*cos(theta);
            radius*sin(theta);
            zeros(1,n_elements)];

    % single target signal
    n = 0:N_snapshots-1;
    s = exp(1j*2*pi*0.05*n);

    X = cell(1,N_rx);
    theta_true = zeros(N_rx,1);

    % simulate received signal at each receiver
    for j = 1:N_rx
        
        % true angle of arrival from receiver j to target
        direction = target - r(:,j);
        theta_true(j) = atan2(direction(2), direction(1));
        
        % unit direction vector
        u = [cos(theta_true(j)); sin(theta_true(j)); 0];
        
        % array manifold vector
        a = exp(-1j* k*(rUCA.'*u));
        
        % received clean signal
        x_clean = a*s;
        
        % add noise
        signal_power = mean(abs(x_clean(:)).^2);
        noise_power = signal_power / (10^(SNR_dB/10));
        noise = sqrt(noise_power/2) * (randn(size(x_clean)) + 1j*randn(size(x_clean)));

        X{j} = x_clean + noise;

    end

    % MUSIC DOA estimation
    phi = linspace(0, 2*pi, 36001);
    theta_hat = NaN(N_rx,1);
    P_music = cell(N_rx,1);
    
    % at each receiver
    for j = 1:N_rx
        
        x = X{j};
        
        % covariance matrix
        Rxx = (x*x') / size(x,2);
        
        % eigendecomposition
        [E,D] = eig(Rxx);
        [~,idx] = sort(real(diag(D)), 'descend');
        E = E(:,idx);
        
        % everything but first eigenvector is noise subspace
        En = E(:, 2:end);
    
        P = zeros(size(phi));
    
        for nidx = 1:length(phi)
    
            u_scan = [cos(phi(nidx)); sin(phi(nidx)); 0];
            a_scan = exp(-1j * k * (rUCA.' * u_scan));

            % MUSIC formula
            P(nidx) = 1 / abs(a_scan' * En * En' * a_scan);
    
        end
    
        P_dB = 10*log10(P / max(P));
        P_music{j} = P_dB;
        
        % choose angle with largest peak for MUSIC spectrum
        [~, loc] = max(P_dB);
        theta_hat(j) = phi(loc);
    
    end

    % Association stage
    rho_hat = estimate_rho(theta_hat, r);

    % Metric fusion stage
    r_hat = metric_fusion(theta_hat, rho_hat, r);

    % print output
    fprintf('True target position: x = %.2f m, y = %.2f m\n', target(1), target(2));
    fprintf('Estimated target position: x = %.2f m, y = %.2f m\n', r_hat(1), r_hat(2));
    fprintf('Position error: %.2f m\n', norm(r_hat - target));

end


function rho_hat = estimate_rho(theta_hat, r)

    % rho12, rho23, rho34, rho41
    pairs = [
        1 2;
        2 3;
        3 4;
        4 1
    ];

    N_rx = size(r,2);
    G = zeros(N_rx,N_rx);
    d = zeros(N_rx,1);

    for q = 1:N_rx

        i = pairs(q,1);
        j = pairs(q,2);

        % distance between two receivers
        baseline_ij = r(:,j) - r(:,i);
        rho_ij = norm(baseline_ij);

        d(q) = rho_ij;

        e_ij = baseline_ij / rho_ij;
        e_ji = -e_ij;

        % unit vectors
        ui = [cos(theta_hat(i)); sin(theta_hat(i))];
        uj = [cos(theta_hat(j)); sin(theta_hat(j))];

        G(q,i) = dot(ui, e_ij);
        G(q,j) = dot(uj, e_ji);

    end

    rho_hat = G \ d;
end


function r_hat = metric_fusion(theta_hat, rho_hat, r)

    N_rx = size(r,2);

    H = kron(ones(N_rx,1), eye(2));
    b = zeros(2*N_rx,1);

    for i = 1:N_rx
        b(2*i-1) = r(1,i) + rho_hat(i)*cos(theta_hat(i));
        b(2*i) = r(2,i) + rho_hat(i)*sin(theta_hat(i));
    end

    r_hat = pinv(H) * b;
end