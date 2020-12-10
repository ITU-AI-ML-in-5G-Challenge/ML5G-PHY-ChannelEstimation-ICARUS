function Psi = dictionary_angles(RXanglesBar,TXanglesBar,RXanglesTilde,TXanglesTilde,Nr,Nt)

Gr = length(RXanglesBar);
Gt = length(TXanglesBar);

AtildeRX = zeros(Nr,Gr);
AtildeTX = zeros(Nt,Gt);

for nr = 1:Nr
    AtildeRX(nr,:) = 1/sqrt(Nr)*exp(1i*pi*(nr-1)*cos(RXanglesBar+RXanglesTilde));
end

for nt = 1:Nt
    AtildeTX(nt,:) = 1/sqrt(Nt)*exp(1i*pi*(nt-1)*cos(TXanglesBar+TXanglesTilde));
end

Psi = kron(conj(AtildeTX),AtildeRX);  % size: (Nt*Nr , Gt*Gr)