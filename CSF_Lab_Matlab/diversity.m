%************ Transmission Diversity and reception Diversity **************
%************      with MPSK or MQAM constellations          **************
%*********** For Matlab 2009, 2010, 2011 and 2012 *************************
frmLen = 100;       % frame length
numPackets = 10000;  % number of packets
EbNo = 0:2:20;      % Eb/No varying to 20 dB
Ntx = 2;              % maximum number of Tx antennas
Nrx = 2;              % maximum number of Rx antennas


% Create a local random stream to be used by random number generators 
hStr = RandStream('mt19937ar', 'Seed', 55408);

% Create MQAM mod-demod objects
%{
P=16;               % modulation order
qammod = modem.qammod('M', 64, 'SymbolOrder', 'Gray');
qamdemod = modem.qamdemod(qammod);
%}


% Create BPSK mod-demod objects
P = 2;				% modulation order
bpskmod = modem.pskmod('M', P, 'SymbolOrder', 'Gray');
bpskdemod = modem.pskdemod(bpskmod);

% Pre-allocate variables for speed
tx2 = zeros(frmLen, Ntx); 
H  = zeros(frmLen, Ntx, Nrx);

r_mrc  = zeros(frmLen, 2);%r_mrc (h1s1+ n1) and (h2s1+n2)
z_mrc = zeros(frmLen, Nrx);%h1*(h1s1+ n1) and h2*(h2s1+n2)
% Alamouti 2tx 1rx
r_alam = zeros(frmLen, 1);% r_alam alamouti received signal 2tx 1 rx 
z_alam2x1 = zeros(frmLen, 1); 
z_1 = zeros(frmLen/Ntx, 1); z_2 = z_1;
%Alamouti 2tx 2rx
r_alam1 = zeros(frmLen, 1);% r_alam alamouti received signal 2tx 1 rx 
r_alam2 = zeros(frmLen, 1);% r_alam alamouti received signal 2tx 1 rx 
z_alam2x2 = zeros(frmLen, 1); 
z_1 = zeros(frmLen/Ntx, 1); z_2 = z_1;z_3 = z_1;z_4 = z_1;


error_wd = zeros(1, numPackets); BERwd = zeros(1, length(EbNo));% without diversity
error_alam2x1= error_wd; BER_alam2x1 = BERwd; error_mrc = error_wd; BER_mrc = BERwd; BERtpsk = BERwd; BERtqam = BERwd;
error_alam2x2= error_wd; BER_alam2x2 = BERwd;
% Set up a figure for visualizing BER results
h = gcf; grid on; hold on;
set(gca, 'yscale', 'log', 'xlim', [EbNo(1), EbNo(end)], 'ylim', [1e-4 1]);
xlabel('Eb/No (dB)'); ylabel('BER'); set(h,'NumberTitle','off');
set(h, 'renderer', 'zbuffer'); set(h,'Name','Transmit vs. Receive Diversity');
title('Transmit vs. Receive Diversity');

% Loop over several EbNo points
for idx = 1:length(EbNo)
    % Loop over the number of packets
    for packetIdx = 1:numPackets
        data = randi(hStr, [0 P-1], frmLen, 1);   % data vector per user and per channel
        tx = modulate(bpskmod, data);             % BPSK modulation it is assumed by default interger for 'inputtype' property
        %tx = modulate(qammod, data);             % QAM modulation it is assumed by default interger for 'inputtype' property
        
        
        % Create the Rayleigh distributed channel response matrix
        %   for two transmit and two receive antennas
        H(1:Ntx:end, :, :) = (randn(hStr, frmLen/2, Ntx, Nrx) + 1i*randn(hStr, frmLen/2, Ntx, Nrx))/sqrt(2);
                  
        %   assume held constant for 2 symbol periods
        H(2:Ntx:end, :, :) = H(1:Ntx:end, :, :); %the channel remains invariant along two consecutive symbols
                                               % H(:,1,1) <->h1 and  H(:,2,1)<-> h2
        
        %*******************************************************************************************
        %******************  Receiver and combiners ************************************************
         %******************************************************************************************
         
       
        
        
        %******************  Without diversity *****************************************************      
        % Received signals
        %   for uncoded 1x1 system
        rwd = awgn(H(:, 1, 1).*tx, EbNo(idx),'measured', hStr);
        r_wd = rwd.*conj(H(:, 1, 1));
        % ML Detector (minimum Euclidean distance)
        demod_wd = demodulate(bpskdemod, r_wd);
        
        
        %********************* MRRC ****************************************************************
        % MRRC Classical Maximum Receive Ratio Combiner 1 x 2
        % calculation of h1s1+ n1 and h2s1+n2
        for i = 1:Nrx
            r_mrc(:, i) = awgn(H(:, 1, i).*tx, EbNo(idx), 0, hStr);
        end
        
        %   Calculation of h1*(h1s1+ n1) and h2*(h2s1+n2)
        for i = 1:Nrx
            z_mrc(:, i) = r_mrc(:, i).* conj(H(:, 1, i));
        end
        % Output at the MRRC combiner
        MRC_comb = sum(z_mrc,2);
        % ML Detector (minimum Euclidean distance)
        demod_mrc = demodulate(bpskdemod, MRC_comb);
        %********************** end of MRRC ********************************************************
        %*******************************************************************************************
        
        
        
        %*******************************************************************************************
        % Alamouti Space-Time Block Encoder, G2, full rate
        %   G2 = [s1 s2; -s2* s1*]
        s1 = tx(1:Ntx:end);
        s2 = tx(2:Ntx:end);
        tx2(1:Ntx:end, :) = [s1 s2];
        tx2(2:Ntx:end, :) = [-conj(s2) conj(s1)]; %generate transmission matrix  [s1 s2; -s2* s1*] 

         %************************* Alamouti 2x1 ***************************************************
        %   for G2-coded 2x1 system - with normalized Tx power, i.e., the
        %	total transmitted power is assumed constant
        r_alam = awgn(sum(H(:, :, 1).*tx2, 2)/sqrt(2), EbNo(idx), 'measured', hStr);
        %pause;
        %  calculates [h1s1+h2s2+n1;-h1s1*+h2s1*+n2] with
        %  r1=h1s1+h2s2+n1 and r2=-h1s1*+h2s1*+n2
        
        % Front-end Combiners - it is assumed that channel response is known at Rx
        hidx = 1:Ntx:length(H);
        z_1 = r_alam(1:Ntx:end).* conj(H(hidx, 1, 1)) +conj(r_alam(2:Ntx:end)).* H(hidx, 2, 1);% h1*r1+h2r2*
        z_2 = r_alam(1:Ntx:end).* conj(H(hidx, 2, 1)) - conj(r_alam(2:Ntx:end)).* H(hidx, 1, 1);%h2*r1-h1r2*
        z_alam2x1(1:Ntx:end) = z_1; z_alam2x1(2:Ntx:end) = z_2;
        % ML Detector (minimum Euclidean distance)
        demod_alam2x1 = demodulate(bpskdemod, z_alam2x1);
        %*******************************************************************************************
        %********************  end Alamouti ********************************************************
       
        %************************** Alamouti 2x2  **************************************************
        %   for G2-coded 2x1 system - with normalized Tx power, i.e., the
        %	total transmitted power is assumed constant
        %  calculates [h1s1+h2s2+n1;-h1s1*+h2s1*+n2] with r1=h1s1+h2s2+n1 and r2=-h1s1*+h2s1*+n2
        r_alam1 = awgn(sum(H(:, :, 1).*tx2, 2)/sqrt(2), EbNo(idx), 'measured', hStr);
        %  calculates [h3s1+h4s2+n3;-h3s2*+h4s1*+n4] with r3=h3s1+h4s2+n3 and r2=-h3s2*+h4s1*+n4
        r_alam2 = awgn(sum(H(:, :, 2).*tx2, 2)/sqrt(2), EbNo(idx), 'measured', hStr);
        %pause
        
        
        % Front-end Combiners - it is assumed that channel response is known at Rx
        hidx = 1:Ntx:length(H);
        z_1 = r_alam1(1:Ntx:end).* conj(H(hidx, 1, 1)) +conj(r_alam1(2:Ntx:end)).* H(hidx, 2, 1);% h1*r1+h2r2*
        z_2 = r_alam1(1:Ntx:end).* conj(H(hidx, 2, 1)) - conj(r_alam1(2:Ntx:end)).* H(hidx, 1, 1);%h2*r1-h1r2*
        z_3 = r_alam2(1:Ntx:end).* conj(H(hidx, 1, 2)) +conj(r_alam2(2:Ntx:end)).* H(hidx, 2, 2);% h3*r3+h4r4*
        z_4 = r_alam2(1:Ntx:end).* conj(H(hidx, 2, 2)) - conj(r_alam2(2:Ntx:end)).* H(hidx, 1, 2);%h4*r3-h3r4*
        z_alam2x2(1:Ntx:end) = z_1+z_3; z_alam2x2(2:Ntx:end) = z_2+z_4;
        % ML Detector (minimum Euclidean distance)
        demod_alam2x2 = demodulate(bpskdemod, z_alam2x2);
        %******************************************************************************************
        %********************  end Alamouti *******************************************************
        

        % Determine errors
        error_wd(packetIdx) = biterr(demod_wd, data);
        error_mrc(packetIdx) = biterr(demod_mrc, data);
        error_alam2x1(packetIdx) = biterr(demod_alam2x1, data);
        error_alam2x2(packetIdx) = biterr(demod_alam2x2, data);
    end % end of FOR loop for numPackets

    % Calculate BER for current idx
    %   for uncoded 1x1 system
    BER_wd(idx) = sum(error_wd)/(numPackets*frmLen);

    %   for G2 coded 2x1 system
    BER_alam2x1(idx) = sum(error_alam2x1)/(numPackets*frmLen);

    %   for G2 coded 2x1 system
    BER_alam2x2(idx) = sum(error_alam2x2)/(numPackets*frmLen);

    
    %   for Maximal-ratio combined 1x2 system
    BER_mrc(idx) = sum(error_mrc)/(numPackets*frmLen);

    %theoretical performance of second-order diversity
    BERtpsk(idx) = berfading(EbNo(idx), 'psk', P, 2);% theoretical BER for M-PSK with fading channel and 2 order diversity 
    
    %BERtqam = berfading(EbNo(idx),'qam',P,2);% theoretical BER for fading M-QAM with fading channel and 2 order diversity 
    
   % %{
    
    % Plot results
    semilogy(EbNo(1:idx), BER_wd(1:idx), 'rs--', ...
             EbNo(1:idx), BER_alam2x1(1:idx), 'ko-', ...
             EbNo(1:idx), BER_alam2x2(1:idx), 'b*--', ...
             EbNo(1:idx), BER_mrc(1:idx), 'gs-.', ...
             EbNo(1:idx), BERtpsk(1:idx), 'mv-');
    legend('No Diversity (1Tx, 1Rx)', 'Alamouti (2Tx, 1Rx)','Alamouti (2Tx, 2Rx)','Maximal-Ratio Combining (1Tx, 2Rx)', ...
           'Theoretical 2nd-Order Diversity');
    set(gca, 'yscale', 'log', 'xlim', [EbNo(1), EbNo(end)], 'ylim', [1e-4 1]);   
    %%}
    %pause;
    drawnow;
    
end  % end of for loop for EbNo
pause;
clf;

semilogy(EbNo, BER_wd,'rs-',EbNo, BER_alam2x1,'ko-', EbNo, BER_alam2x2,'b*--', EbNo, BER_mrc,'gs-.');
legend('No Diversity (1Tx, 1Rx)', 'Alamouti (2Tx, 1Rx)','Alamouti (2Tx, 2Rx)','Maximal-Ratio Combining (1Tx, 2Rx)', ...
           'Theoretical 2nd-Order Diversity');
set(gca, 'yscale', 'log', 'xlim', [EbNo(1), EbNo(end)], 'ylim', [1e-4 1]);
pause;

clf
