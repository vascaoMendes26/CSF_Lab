clear all;
close all; 
clc;
mode = '256QAM'; % 'QPSK', '64QAM' ou '256QAM'
CHANNEL='AWGN'; % 'AWGN', 'RAYL', 'REAL', 'XTAP', 'RRND'
SYSTEM = 'OFDM'; % 'OFDM', 'SCFD'
Ssat_value = 5;
L=4; % L-th order diversity
FASE = 2; %Fase 1 s/Amplificador Fase 2 c/Amplificador
p_rapp = 3; %p_rapp = 100->Ideal p_rapp = 1->N_Ideal P=2 p=10 e p=50

if strcmp(mode,'QPSK')
    M      = 4; % 4QPSK
    Eb     = 1;
elseif strcmp(mode,'64QAM')
    M      = 64; %64QAM
    Eb     = 7;
elseif strcmp(mode, '256QAM')
    M      = 256; %256QAM
    Eb     = 21.25;
end
   
Nbs    = log2(M); % Calcula o Num de bits por simbolo, importante para calcular a potencia 

EN     = [-5:2:50]'+0*100; en = 10 .^(EN/10) ; % Aqui podemos mudar o valor de SNR que queremos ter mudando o array

N=512;
NSlot=1000; %Variar para gráfico mais curvado, aumentamos numero iteracoes
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
    
    % --- GERAÇÃO DO CANAL ---
    if strcmp(CHANNEL,'REAL')
        Hk=zeros(N,L); 
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
            end;
        end; 
    elseif strcmp(CHANNEL,'XTAP')
        Hk=zeros(N,L); 
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
            end;
        end;
    elseif strcmp(CHANNEL,'RRND')
        Hk=zeros(N,L); 
        tau=rand(NRay,1)*Tg;
            for l=1:L
                alpha=ones(NRay,1).*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2*NRay);
                for nRay=1:NRay
                    Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
                end;
            end; 
    elseif strcmp(CHANNEL,'RAYL')
        Hk=(randn(N,L)+j*randn(N,L))/sqrt(2);
    elseif strcmp(CHANNEL,'AWGN')
        Hk=ones(N,L).*exp(j*2*pi*rand(N,L));
    end;
    H2k=abs(Hk).^2;
    if (L==1) sH2k=H2k; else sH2k=sum(H2k')'; end;
    
    % --- TRANSMISSOR ---
    Tx_Data = randi([0 M-1], N, 1);      % Gera inteiros
    Ak_Tx = qammod(Tx_Data, M);          % Modula corretamente
    
    % --------------------------------------------------

    if (FASE == 1)
        % --- FASE 1: AMPLIFICADOR IDEAL ---
        if strcmp(SYSTEM, 'OFDM')
            an_Tx = fftshift(ifft(fftshift(Ak_Tx)));
            Tx_To_Channel = fftshift(fft(fftshift(an_Tx))); 
        elseif strcmp(SYSTEM, 'SCFD')
            an_Tx = Ak_Tx;
            Tx_To_Channel = fftshift(fft(fftshift(an_Tx))); 
        end 

    elseif (FASE == 2)
        % --- FASE 2: AMPLIFICADOR REAL (RAPP) ---
        if strcmp(SYSTEM, 'OFDM') 
            
            an_Tx=fftshift(ifft(fftshift(Ak_Tx))); 
            
            % Amplificador é no tempo
            amplitude_tx = abs(an_Tx);
            phase_tx = angle(an_Tx);

            Ssat = Ssat_value; 
            amplitude_tx_mean = mean(amplitude_tx); 
            p = p_rapp; 
            satlevel = amplitude_tx_mean*10^(Ssat/10);

            A = amplitude_tx./(1+(amplitude_tx./satlevel).^(2*p)).^(1/(2*p));
            an_Tx_Distorted = A.*exp(j*phase_tx); 
            
            % Volta para Frequência 
            Tx_To_Channel = fftshift(fft(fftshift(an_Tx_Distorted)));

        elseif strcmp(SYSTEM, 'SCFD') 
            
            % SC-FDE: Símbolos já estão no Tempo
            an_Tx = Ak_Tx;

            amplitude_tx = abs(an_Tx);   
            phase_tx = angle(an_Tx);     

            Ssat = Ssat_value; 
            amplitude_tx_mean = mean(amplitude_tx); 
            p = p_rapp; 
            satlevel = amplitude_tx_mean*10^(Ssat/10);

            A = amplitude_tx./(1+(amplitude_tx./satlevel).^(2*p)).^(1/(2*p));
            an_Tx_Distorted = A.*exp(j*phase_tx); 
            
            % Volta para Frequência 
            Tx_To_Channel = fftshift(fft(fftshift(an_Tx_Distorted)));
            
        else
            error(['ERRO CRÍTICO: O sistema "', SYSTEM, '" não existe.']);
        end
    end
    
    % --- CANAL E RECETOR ---
     for nEN=1:NEN
        Yk=zeros(N,L);
        for l=1:L
            Yk(:,l)=Tx_To_Channel.*Hk(:,l)+(randn(N,1)+j*randn(N,1))*sigma(nEN);
        end;
        YIk=0;
        for l=1:L
            YIk = YIk +Yk(:,l).*conj(Hk(:,l));
        end;
        YIk=YIk./sH2k; % Sinal Equalizado na Frequência
        
        % --- BLOCO DO RECETOR ---
        if strcmp(SYSTEM, 'OFDM')
            Rx_Symbols = YIk;
        elseif strcmp(SYSTEM, 'SCFD')
            Rx_Symbols = fftshift(ifft(fftshift(YIk)));
        end
        
        Rx_Data = qamdemod(Rx_Symbols, M);                   
        [numBitsErr, ~] = biterr(Tx_Data, Rx_Data, log2(M)); 
        NErr(nEN,1)=NErr(nEN,1)+numBitsErr;
    end;

    if (rem(nn,100)==0) nn, end;
end;

% BER in Rayleigh channel and L-branch diversity 
aux=sqrt(en./(1+en));Pb_tr=0;
for l=0:L-1
    Pb_tr=Pb_tr+combin(L-1+l,l)*((1+aux)/2).^l;
end;
Pb_tr=Pb_tr.*((1-aux)/2).^L;

% BER in AWGN channel
PbAWGN=q_x(sqrt(2*L*en));

Pb=NErr/NSlot/N/Nbs;

%figure;
semilogy(EN,Pb,'g-*',EN,PbAWGN,'b:',EN,Pb_tr,'b*:')
xlabel('E_b/N_0(dB)'),ylabel('BER')
axis([-5 50 1e-4 1])
%pause,clf;

%--------------
%Explicação Plots: 
%Linha (Azul fina): "O ideal modulação simples (QPSK)."
%Linha (Verde): "O sistema MQAM a funcionar. Precisa de mais sinal que o QPSK (Azul), mas desce rápido se o canal é estável (AWGN)."
%Linha da Direita (Azul estrelada):Rayleigh.  mostrar como o canal seria mau se não fosse AWGN."
%---------------
%Explicação de Variaveis:
%Ak -> Sempre que temos Ak tamos a passar amplittude para freq;
%An -> Sempre que temos An tamos a passar amplitutde para tempo;
%---------------