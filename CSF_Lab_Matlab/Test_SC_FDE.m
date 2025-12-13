% Simulates SC-FDE with IB-DFE
EN=[-6:2:16]'+0*100; 

en = 10 .^(EN/10) ;

N=256;
CHANNEL='XTAP';
NSlot=100;
L=1; % L-th order diversity
NIter=5;
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

Eb=N;

sigma=sqrt(Eb/2 ./en);
NSR=1/2 ./(en); 

NEN=length(EN);
NErr=zeros(NEN,NIter);
PbMFB=zeros(NEN,1);
Rho_Tot=zeros(NEN,NIter);

for nn=1:NSlot

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
    % Matched filter bound for that channel
    PbMFB=PbMFB+q_x(sqrt(2*en*mean(sH2k)));
    
    an_Tx=sign(randn(N,1))+j*sign(randn(N,1));
    Ak_Tx=fftshift(fft(fftshift(an_Tx)));
    
    for nEN=1:NEN
        
        Yk=zeros(N,L);
        for l=1:L
            Yk(:,l)=Ak_Tx.*Hk(:,l)+(randn(N,1)+j*randn(N,1))*sigma(nEN);
        end;
        
        Rho=0; Ak_est=zeros(N,1);
        Fk=zeros(N,L);
        for nIter=1:NIter
            
            % Fk
%            Ik= (1-Rho^2)*sH2k + NSR(nEN) ;
            Ik= (1-Rho^2)*sH2k + NSR(nEN) ;
            for l=1:L
                Fk(:,l) = conj(Hk(:,l))./Ik;
            end;
            
            % gamma
            gamma=0;
            for l=1:L
                gamma=gamma+sum(Fk(:,l).*Hk(:,l));
            end;
            gamma=gamma/N;
            
            % Normalization
            Fk=Fk/gamma;
            gamma=1;
            
            % Bk
            Bk=0;
            for l=1:L
                Bk=Bk+Fk(:,l).*Hk(:,l);
            end;
            Bk=Rho*(Bk-gamma);
            
            % IB-DFE
            YIk=0;
            for l=1:L
                YIk = YIk +Yk(:,l).*Fk(:,l) ;
            end;
            YIk=YIk-Bk.*Ak_est;
            yIn = fftshift(ifft(fftshift(YIk))) ;
        
            an_est=sign(real(yIn))+j*sign(imag(yIn));
            Ak_est=fftshift(fft(fftshift(an_est)));
            aux = sum( abs(sign(real(an_Tx)-real(an_est))) + ...
                    abs(sign(imag(an_Tx)-imag(an_est))) ) ;
            NErr(nEN,nIter)=NErr(nEN,nIter)+aux;
            
            Rho=sum(an_est.*conj(an_Tx))/sum(abs(an_Tx).^2);
            Rho=abs(Rho);
            Rho_Tot(nEN,nIter)=Rho_Tot(nEN,nIter)+Rho;
        end;
        
    end;

    if (rem(nn,100)==0)
        nn
    end;
end;

% BER in Rayleigh channel and L-branch diversity [Proakis]
aux=sqrt(en./(1+en));
Pb_tr=0;
L=1;
for l=0:L-1
    Pb_tr=Pb_tr+combin(L-1+l,l)*((1+aux)/2).^l;
end;
Pb_tr=Pb_tr.*((1-aux)/2).^L;

% BER in AWGN channel
Pb_ta=q_x(sqrt(2*en));

Pb=NErr/2/N/NSlot;
Rho_Tot=Rho_Tot/NSlot;
PbMFB=PbMFB/NSlot;

plot(EN,Rho_Tot,'-')
xlabel('E[E_b]/N_0(dB)'),ylabel('\rho')
axis([0 20 0 1])
pause,clf;

%figure;
semilogy(EN,Pb,'-',EN,PbMFB,':')
xlabel('E_b/N_0(dB)'),ylabel('BER')
axis([0 20 1e-4 1])
%pause,clf;
