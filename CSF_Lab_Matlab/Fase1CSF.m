mode = '64QAM'; % Choose modulation mode
CHANNEL='AWGN';%AWGN or RAYL 
L=1; % L-th order diversity of antennas
Amp = 1;

if strcmp(mode,'QPSK')
    M      = 4; % 4 QPSK
    Eb     = 1;
    
elseif strcmp(mode,'64QAM')
    M      = 64; %64QAM
    Eb     = 7;

elseif strcmp(mode, '256QAM')
    M      = 256; %256QAM
    Eb     = 21.25;
end
   
Nbs    = log2(M); % Num de bits por simbolo

EN     = [-5:2:50]'+0*100; en = 10 .^(EN/10) ; % SNR 

N=512; %FFT Transform 
NSlot=1000; % N samples of graph (Monte Carlo)
Ts=4e-6; % Block duration
Tg=0.2*Ts; % Cyclic prefix durration
f=[-N/2:N/2-1]'/Ts; % frequencies

if (CHANNEL=='AWGN')
    alfa_med=1;tau=0;NRay=1;
elseif (CHANNEL=='REAL')
%    load pdp.940
%    load PDP_ChA.dat,pdp=PDP_ChA;
    load pdp_hipc.dat,pdp=pdp_hipc;
    tau=pdp(:,3);
    alpha_med=10 .^((pdp(:,1))/20);
    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));
    NRay=length(alpha_med);
%    NRay=10;tau=[0:1:NRay-1]'*0.1/NRay*Ts; % Tg=0.1Ts
elseif (CHANNEL=='XTAP')
    NRay=32;
    tau=[0:NRay-1]'*Ts/N;
    alpha_med=ones(NRay,1);
    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));
end;

sigma=sqrt(Eb/2 ./en); 
NSR=1/2 ./(en); 
NEN=length(EN);
NErr=zeros(NEN,1);


for nn=1:NSlot
    
    %rand('state',nn*1234567); randn('state',nn*1234567);
    % This means the same channel for each slot

    if (CHANNEL=='REAL')
        Hk=zeros(N,L); 
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
            end;
        end; 
    elseif (CHANNEL=='XTAP')
        Hk=zeros(N,L); 
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
            end;
        end;
    elseif (CHANNEL=='RRND')
        Hk=zeros(N,L); 
        tau=rand(NRay,1)*Tg;
            for l=1:L
                alpha=ones(NRay,1).*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2*NRay);
                for nRay=1:NRay
                    Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
                end;
            end; 
    elseif (CHANNEL=='RAYL')
        Hk=(randn(N,L)+j*randn(N,L))/sqrt(2);
    elseif (CHANNEL=='AWGN')
        Hk=ones(N,L).*exp(j*2*pi*rand(N,L));
    end;
    H2k=abs(Hk).^2;
    if (L==1) sH2k=H2k; else sH2k=sum(H2k')'; end;
    

    k_root = sqrt(M); 
    % Gera niveis PAM: -7, -5, ... 5, 7 (exemplo para 64QAM)
    Re = 2*floor(k_root*rand(N,1)) - (k_root-1);
    Im = 2*floor(k_root*rand(N,1)) - (k_root-1);
    Ak_Tx_Raw = Re + j*Im;
    
    % Normalizar para potência unitária (importante para o amplificador funcionar bem)
    scale_factor = sqrt(mean(abs(Ak_Tx_Raw).^2)); 
    Ak_Tx = Ak_Tx_Raw / scale_factor; 

    % Passar para o Tempo (OFDM)
    an_Tx = fftshift(ifft(fftshift(Ak_Tx)));
    
    %Cálculo das Amplificações
    if(Amp = 1) %Para não separar código aqui definimos se fazemos fase 2 ou não
        P_SSPA = 2;
        V_sat  = 1.0;
        

    for nEN=1:NEN
        Yk=zeros(N,L);
        for l=1:L
            Yk(:,l)=Ak_Tx.*Hk(:,l)+(randn(N,1)+j*randn(N,1))*sigma(nEN);
        end;
        YIk=0;
        for l=1:L
            YIk = YIk +Yk(:,l).*conj(Hk(:,l));
        end;
        YIk=YIk./sH2k;

        Ak_Rx=sign(real(YIk))+j*sign(imag(YIk));
        aux = sum( abs(real(Ak_Tx)-real(Ak_Rx)) + ...
                abs(imag(Ak_Tx)-imag(Ak_Rx)) ) / 2 ;
        NErr(nEN,1)=NErr(nEN,1)+aux;
    end;

    if (rem(nn,100)==0) nn, end;
end;

% BER in Rayleigh channel and L-branch diversity [Proakis]
aux=sqrt(en./(1+en));Pb_tr=0;
for l=0:L-1
    Pb_tr=Pb_tr+combin(L-1+l,l)*((1+aux)/2).^l;
end;
Pb_tr=Pb_tr.*((1-aux)/2).^L;

% BER in AWGN channel
PbAWGN=q_x(sqrt(2*L*en));

Pb=NErr/NSlot/N/2;

%figure;
semilogy(EN,Pb,'g-*',EN,PbAWGN,'b:',EN,Pb_tr,'b*:')
xlabel('E_b/N_0(dB)'),ylabel('BER')
axis([-5 50 1e-4 1])
%pause,clf;