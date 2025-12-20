clear all;
close all; 
clc;

% --- CONFIGURAÇÕES ---
mode = '256QAM';    % 'QPSK', '64QAM' ou '256QAM'
CHANNEL = 'RAYL';   % 'AWGN', 'RAYL', 'REAL', 'XTAP', 'RRND'
SYSTEM = 'OFDM';    % 'OFDM', 'SCFD'
Ssat_value = 5;     % Nível de Saturação (quanto menor, mais distorção)
L = 1;              % Diversidade na receção
FASE = 2;           % 1 = Ideal, 2 = Com Amplificador (Rapp)
p_rapp = 1;         % 1 = Forte Distorção, 10+ = Linear
IBO_dB = 6;         % Input Backoff: Mete 0, 3, 6, 10 para testar

% --- CONFIGURAÇÃO DA MODULAÇÃO ---
if strcmp(mode,'QPSK')
    M = 4; Eb = 1;
elseif strcmp(mode,'64QAM')
    M = 64; Eb = 7;
elseif strcmp(mode, '256QAM')
    M = 256; Eb = 21.25;
end
   
Nbs = log2(M); 
EN = [-5:2:70]'; 
en = 10.^(EN/10); 

N = 512;
NSlot = 1000;       % Numero de iterações
Ts = 4e-6; 
Tg = 0.2*Ts; 
f = [-N/2:N/2-1]'/Ts;

% --- CONFIGURAÇÃO DO CANAL ---
if (CHANNEL=='AWGN')
    alpha_med=1; tau=0; NRay=1;
elseif (CHANNEL=='REAL')
    load pdp_hipc.dat, pdp=pdp_hipc;
    tau=pdp(:,3);
    alpha_med=10.^((pdp(:,1))/20);
    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));
    NRay=length(alpha_med);
elseif (CHANNEL=='XTAP')
    NRay=32;
    tau=[0:NRay-1]'*Ts/N;
    alpha_med=ones(NRay,1);
    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));
end

sigma = sqrt(Eb/2 ./en); 
NEN = length(EN);
NErr = zeros(NEN,1);
scale_factor = 1; % Inicialização para não dar erro na FASE 1

% --- LOOP PRINCIPAL ---
for nn=1:NSlot
    
    % --- GERAÇÃO DO CANAL ---
    if strcmp(CHANNEL,'REAL') || strcmp(CHANNEL,'XTAP')
        Hk=zeros(N,L); 
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+1j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-1j*2*pi*f*tau(nRay));
            end
        end 
    elseif strcmp(CHANNEL,'RAYL')
        Hk=(randn(N,L)+1j*randn(N,L))/sqrt(2);
    elseif strcmp(CHANNEL,'AWGN')
        Hk=ones(N,L).*exp(1j*2*pi*rand(N,L));
    end
    
    H2k=abs(Hk).^2;
    if (L==1) sH2k=H2k; else sH2k=sum(H2k')'; end
    
    % --- TRANSMISSOR (DADOS) ---
    Tx_Data = randi([0 M-1], N, 1);
    Ak_Tx = qammod(Tx_Data, M); 
    
    % --- PROCESSAMENTO DO SINAL (FASE 1 vs FASE 2) ---
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
        % --- FASE 2: AMPLIFICADOR REAL (RAPP) COM BACKOFF ---
        
        % 1. Obter o sinal no Domínio do Tempo (Original)
        if strcmp(SYSTEM, 'OFDM')
            an_Tx_Original = fftshift(ifft(fftshift(Ak_Tx))); 
        elseif strcmp(SYSTEM, 'SCFD')
            an_Tx_Original = Ak_Tx; 
        end

        % 2. Calcular Nível de Saturação (Baseado no sinal SEM Backoff)
        amp_ref = abs(an_Tx_Original);
        fixed_sat_level = mean(amp_ref) * 10^(Ssat_value/10);

        % 3. Aplicar Input Backoff (Reduzir amplitude entrada)
        scale_factor = 10^(-IBO_dB/20); 
        an_Tx_Input = an_Tx_Original * scale_factor;

        % 4. Modelo de Rapp (Distorção)
        amplitude_tx = abs(an_Tx_Input);
        phase_tx = angle(an_Tx_Input);
        
        % Fórmula Rapp 
        A = amplitude_tx ./ (1 + (amplitude_tx./fixed_sat_level).^(2*p_rapp)).^(1/(2*p_rapp));
        an_Tx_Distorted = A .* exp(1j*phase_tx);
        
        % 5. Voltar para Frequência para o Canal
        Tx_To_Channel = fftshift(fft(fftshift(an_Tx_Distorted)));
    end
    
 
    % --- CANAL E RECETOR ---
     for nEN=1:NEN
        Yk=zeros(N,L);
        for l=1:L
            % Adiciona Ruído
            noise = (randn(N,1)+1j*randn(N,1))*sigma(nEN);
            Yk(:,l) = Tx_To_Channel.*Hk(:,l) + noise;
        end
        
        % Equalização MRC
        YIk=0;
        for l=1:L
            YIk = YIk + Yk(:,l).*conj(Hk(:,l));
        end
        YIk = YIk./sH2k; 
        
        % --- NORMALIZAÇÃO (AGC) ---
        YIk = YIk ./ scale_factor;
        % --------------------------

        % Demodulação
        if strcmp(SYSTEM, 'OFDM')
            Rx_Symbols = YIk;
        elseif strcmp(SYSTEM, 'SCFD')
            Rx_Symbols = fftshift(ifft(fftshift(YIk)));
        end
        
        Rx_Data = qamdemod(Rx_Symbols, M);                    
        [numBitsErr, ~] = biterr(Tx_Data, Rx_Data, log2(M)); 
        NErr(nEN,1) = NErr(nEN,1) + numBitsErr;
    end

    if (rem(nn,100)==0) disp(['Slot: ' num2str(nn)]); end
end

% --- PLOT (CÓDIGO ORIGINAL) ---
% BER in Rayleigh channel and L-branch diversity 
aux=sqrt(en./(1+en)); Pb_tr=0;
for l=0:L-1
    Pb_tr=Pb_tr+combin(L-1+l,l)*((1+aux)/2).^l;
end
Pb_tr=Pb_tr.*((1-aux)/2).^L;

% BER in AWGN channel
PbAWGN=q_x(sqrt(2*L*en));

Pb=NErr/NSlot/N/Nbs;

figure;
semilogy(EN,Pb,'g-*',EN,PbAWGN,'b:',EN,Pb_tr,'b*:')
xlabel('E_b/N_0(dB)'),ylabel('BER')
axis([-5 50 1e-4 1])
grid on;