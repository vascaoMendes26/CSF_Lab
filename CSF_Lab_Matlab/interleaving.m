
%Example: Interleave and deinterleave data with random interleaver
    permVec = randperm(7)';  % Random permutation vector
    hInt = comm.gpu.BlockInterleaver(permVec);
    hDeInt = comm.gpu.BlockDeinterleaver(permVec);
    data = randi(9, 7, 1);
    intData = step(hInt, data);
    deIntData = step(hDeInt, intData);
    % compare the original sequence, interleaved sequence,
    % and restored sequence
    [data, intData, deIntData]