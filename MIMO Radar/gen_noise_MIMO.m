% Generate noise
function [noise, noise_power] = gen_noise_MIMO(signal, SNR_dB)

    % estimate signal power
    signal_power = mean(abs(signal(:)).^2);

    % compute required noise power for desired SNR
    noise_power = signal_power / (10^(SNR_dB/10));

    % generate complex noise
    noise = sqrt(noise_power/2) * ...
        (randn(size(signal)) + 1j * randn(size(signal)));

end