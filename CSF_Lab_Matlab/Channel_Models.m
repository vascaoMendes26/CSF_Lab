%PDP='AWGN'; AWGN channel
%PDP='RAYL'; Uncorrelated fading on the different subcarriers
%PDP='REAL'; Specific Power delay profile (PDP), as in file 
%                               PDP_ChA.dat: Channel Vehicular A from UMTS
%                               PDP_HipC.dat: Channel C from HIPERLAN/2
%                               PDP.940: Typical large indoor scenario
%                               (power of the ray in column 1 and delay in
%                               column 3; the others are not used (now))
%PDP='XTAP'; NRay taps with delays [0:NRay-1]*Ts/N
%PDP='RRND'; NRay taps with random dealys, uniformly distributted in [0:Tg-Ts/N]

Ts=4e-6; % Block duration
Tg=0.2*Ts; % Cyclic prefix durration
f=[-N/2:N/2-1]'*Ts; % frequencies

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
    NRay=6;
    tau=[0:NRay-1]'*Ts/N;
    alpha_med=[0:NRay-1]';
end;    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));


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

    PbMFB=PbMFB+q_x(sqrt(2*en*sum(sum(abs(Hk).^2))/N));
    if (L==1) sH2k=H2k; else sH2k=sum(H2k')'; end;

