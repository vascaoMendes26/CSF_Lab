function bits_hat = DeOFDM(yk_hat,mode)

    % Determine modulation order
    if strcmp(mode, '64QAM')
        M = 64;
    elseif strcmp(mode, 'QPSK')
        M = 4;
    elseif strcmp(mode, '256QAM')
        M = 256;
    else
        error('OFDM: Mode not implemented');
    end

    bits_hat = qamdemod(yk_hat, M, 'OutputType', 'bit');
end

