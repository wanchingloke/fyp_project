% Generate transmit baseband signal for one dwell using m-sequences
function [baseband_signal, mseq_codes] = gen_Tx_signal_MIMO(A, num_PRIs, M_range_bins, Ntx)
    
    m = 7;
    Nc = 2^m - 1;

    all_taps = {
        [7 1]
        [7 3]
        [7 3 2 1]
        [7 4 3 2]
        [7 6 5 4]
        [7 6 3 1]
        [7 5 4 3]
        [7 6 4 2]
    };

    tap_list = all_taps(1:Ntx);     % select first Ntx taps
    mseq_codes = zeros(Ntx, Nc);    % generate Tx codes
    
    for tx = 1:Ntx
        pn = generate_mseq(tap_list{tx}, m);
        code = 2*pn - 1;
        mseq_codes(tx,:) = code(1:Nc);
    end

    total_samples = M_range_bins * Nc * num_PRIs;
    baseband_signal = zeros(Ntx, total_samples);

    for tx = 1:Ntx
        one_pulse = A * mseq_codes(tx,:);

        one_PRI = zeros(1, M_range_bins * Nc);
        one_PRI(1:Nc) = one_pulse;

        baseband_signal(tx,:) = repmat(one_PRI, 1, num_PRIs);
    end
end

function seq = generate_mseq(taps, m)
    reg = ones(1,m);            % initialise shift registers with all ones
    seq = zeros(1, 2^m - 1);    % initialise sequence with all zeros

    % Loop through each chip
    for n = 1:(2^m - 1)
        seq(n) = reg(end);                  % take the output bit from last register 
        feedback = mod(sum(reg(taps)), 2);  % generate feedback bit using XOR
        reg(2:end) = reg(1:end-1);          % shift the register
        reg(1) = feedback;                  % insert new feedback bit at the start   
    end
end