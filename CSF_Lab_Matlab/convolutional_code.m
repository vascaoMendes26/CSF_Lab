cl=3;%constrain lenght
tlenght=12;
codegen=[7 5]; % Polynomials generators
coderate=1/2;
trellis=poly2trellis(cl,codegen);
hstr=RandStream('mt19937ar','seed',65321);
hqpskmod = comm.QPSKModulator('BitInput',true);
cont=1;hqpskdemod = comm.QPSKDemodulator('BitOutput',true);

error=zeros(size(11,1)); Pb=zeros(size(11,1));
for snrdb=1:1:6
    erro=0;nbitstot=0;
    while erro<500
        
        hchannel = comm.AWGNChannel('NoiseMethod',...
    'Signal to noise ratio (SNR)','SNR',snrdb);    %channel noise object (white gaussian noise)

        bits_inf=randi(hstr,[0 1],256,1); %generates the information bits to be sent to the channel
        bits_iterm=[bits_inf(:,:);zeros(cl,1)];% appends the termination bits to the information bits (it must be equal to the constrain length)
        seq_lenght=tlenght+cl;


        bits_code=convenc(bits_iterm,trellis); %encoding of the information bits
        bits_mod=step(hqpskmod,bits_code);  %the coded bits are modulated using an QPSK constellation
        k=0;  
        houtputchannel=step(hchannel,bits_mod); 
        bits_demod=step(hqpskdemod,houtputchannel); %signal demodulation
        bits_decod=vitdec(bits_demod,trellis,seq_lenght,'term','hard'); %viterbi algorithm to decode
        nbitstot=nbitstot+259;
        erro=erro+abs(sum(bits_decod - bits_iterm));
    end
error(cont,1)=erro;    %error
if cont==1
    fprintf('BER with convolutional code in a AWGN channel');
end;
Pb(cont,1)=(erro/nbitstot)
cont=cont+1;

end
snr=[1 2 3 4 5 6];
semilogy(snr,Pb);
xlabel('E_b/N_0(dB)'),ylabel('BER')
pause;
clf;