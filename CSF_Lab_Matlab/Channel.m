function yk =Channel(sk,Hk,Mr,Nc);
    
    %xk     -> symbols size (Nc*Mt)
    %Mt     -> Tx antennas 
    %Mr     -> Rx antennas
    %Nc     -> Number of channels
    %yk     -> Symbols out (Nc,Mr)

    yk = zeros(Nc,Mr);
    for n=1:Nc
        yk(n,:) = Hk(:,:,n)*(sk(n,:).');
    end

end

