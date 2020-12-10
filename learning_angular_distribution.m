%Learning the angle distribution

%Size of Hk: (64,16,256)
%Size of Ftr: (16, 200)
%Size of Wtr: (64, 400)
%Size of measurement matrix Phi: (Ntrain*Lr,Nt*Nr)

sett = 1;

load(strcat('training_channels_',num2str(sett),'.mat'));

tugfe = channels_saved;
tugfe2 = sqrt(256)*ifft(tugfe,[],3);
Nc = 1000;
K_total = zeros(Nc,1);
for nc = 1:Nc
    vektor1 = vec(sum(abs(tugfe2(:,:,:,nc)).^2,[1, 2]));
    for k=256:-1:1
        if vektor1(k)>10^(-15)
            K_total(nc) = k;
            break
        end
    end
end

Kmax = max(K_total);


Nt = 16; % Number of TX antennas
Nr = 64; % Number of RX antennas
Lt = 2;  % Number of TX RF chains
Lr = 4;  % Number of RX RF chains
Ns = 2;  % Number of data streams to be transmitted
M = 100; % Number of training symbols to be received for each one of the available channels

Gr = 96;  % AoA grid size
Gt = 24;  % AoD grid size

muuTot = zeros(Gr,Gt,Kmax,Nc);

betaR = 0.01;
betaT = 0.01;

RXanglesBar = linspace(0,pi,Gr);
TXanglesBar = linspace(0,pi,Gt);
RXanglesTilde = zeros(1,Gr);
TXanglesTilde = zeros(1,Gt);
gammaa = 1/0.0001;

Psi = dictionary_angles(RXanglesBar,TXanglesBar,RXanglesTilde,TXanglesTilde,Nr,Nt);
Psi2 = gammaa*(Psi'*Psi);


for nc = 1:Nc
    
    K = K_total(nc);
    
    y = zeros(Nr*Nt,K);
    for k = 1:K
        y(:,k) = vec(tugfe2(:,:,k,nc))+sqrt(0.5*0.0001)*(randn(Nr*Nt,1)+1i*randn(Nr*Nt,1));
    end
    
    
    yy2 = gammaa*Psi'*y;
    
    alpha = ones(Gr,Gt);
    xOld = zeros(Gr,Gt,K);
    diffEM = 100;
    epsilonEM = 10^(-4);
    maxIterEM = 50;
    iterEM = 0;
    
    %Start iterations
    muu = zeros(Gr,Gt,K);
    chii  = repmat(1./alpha, [1,1,K]);
    
    while (iterEM < maxIterEM) && (diffEM > epsilonEM)
        
        
        alphaLeft = zeros(Gr+2,Gt+2);
        alphaRight = zeros(Gr+2,Gt+2);
        alphaLower = zeros(Gr+2,Gt+2);
        alphaUpper = zeros(Gr+2,Gt+2);
        
        alphaLeft(2:end-1,1:end-2) = alpha;
        alphaRight(2:end-1,3:end) = alpha;
        alphaLower(3:end,2:end-1) = alpha;
        alphaUpper(1:end-2,2:end-1) = alpha;
        
        etaa = alpha + betaR*alphaLower(2:end-1,2:end-1) + betaR*alphaUpper(2:end-1,2:end-1)...
            + betaT*alphaLeft(2:end-1,2:end-1) + betaT*alphaRight(2:end-1,2:end-1);
        
        
        temporInv = Psi2+diag(vec(etaa));
        chii = repmat(reshape(diag(temporInv\eye(Gr*Gt)),Gr,Gt),[1,1,K]);
        muu=reshape(temporInv\yy2,Gr,Gt,K);
        %M-Step: Hyper-parameter Learning
        
        %Step 1: Update for for \alpha
        muchiLeft = zeros(Gr+2,Gt+2);
        muchiRight = zeros(Gr+2,Gt+2);
        muchiLower = zeros(Gr+2,Gt+2);
        muchiUpper = zeros(Gr+2,Gt+2);
        
        temporrMatr = sum(abs(muu).^2 + chii,3);
        muchiLeft(2:end-1,1:end-2) = temporrMatr;
        muchiRight(2:end-1,3:end) = temporrMatr;
        muchiLower(3:end,2:end-1) = temporrMatr;
        muchiUpper(1:end-2,2:end-1) = temporrMatr;
        
        omegaa = temporrMatr + betaR*muchiLower(2:end-1,2:end-1) + betaR*muchiUpper(2:end-1,2:end-1)...
            + betaT*muchiLeft(2:end-1,2:end-1) + betaT*muchiRight(2:end-1,2:end-1);
        
        alpha = K./(omegaa);
        
        diffEM = sum(abs(xOld-muu).^2,'all')/sum(abs(muu).^2,'all');
        xOld = muu;
        
        iterEM  = iterEM+1;
        
    end
    nc
    
    sum(abs(fft(reshape(Psi*reshape(muu,96*24,K),64,16,K),256,3)/sqrt(256)-tugfe(:,:,:,nc)).^2,'all')...
        /sum(abs(tugfe(:,:,:,nc)).^2,'all')
    
    muuTot(:,:,1:K,nc) = muu;
    
end

save(strcat('muu_training_',num2str(sett),'.mat'),'muuTot','-v7.3');
