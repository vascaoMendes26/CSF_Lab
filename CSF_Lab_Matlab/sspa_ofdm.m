

    % Simulates OFDM

    N=512;
    NSlot=1000;




    
   
    Ak_Tx=sign(randn(N,1))+j*sign(randn(N,1)); %qpsk 1/-1 + j/-j      
  
    
    An_Tx = fftshift(ifft(fftshift(Ak_Tx))); %ifft
    %conversão para polares
    envelope_tx = abs(An_Tx);   
    phase_tx = angle(An_Tx);
    
   %SSPA
    Ssat=5; %exemplo de valor
    envelope_tx_mean=mean(envelope_tx);
    p = 10;% P pode ser outro P=1, 2, 5, 10 ou 100
    satlevel = envelope_tx_mean*10^(Ssat/10);   
    A = envelope_tx./(1+(envelope_tx./satlevel).^(2*p)).^(1/(2*p));
    Fi = phase_tx;      
    
  
    An_Tx = A.*exp(j*Fi);
    Ak_Tx=fftshift(fft(fftshift(An_Tx)));
   
    
