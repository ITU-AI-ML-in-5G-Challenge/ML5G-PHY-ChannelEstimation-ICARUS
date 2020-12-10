load('learned_parameters.mat')
SNR = 1;  % 1, 2, or 3
M = 20; % 20, 40, or 80, Number of training symbols

if SNR == 1
    gammaLow = 10^(-20/10);
    gammaUpp = 10^(-11/10);
elseif SNR == 2
    gammaLow = 10^(-11/10);
    gammaUpp = 10^(-6/10);
else
    gammaLow = 10^(-6/10);
    gammaUpp = 10^(-0/10);
end




Nt = 16; % Number of TX antennas
Nr = 64; % Number of RX antennas
Lt = 2;  % Number of TX RF chains
Lr = 4;  % Number of RX RF chains

Gr = 96;  % AoA grid size
Gt = 32;  % AoD grid size
betaa = 0.001;

estimated_channels = zeros(1000,Nr,Nt,256);
received_real = h5read(strcat('test_dataset_v3_',num2str(M),'_pilots_',num2str(SNR),'_data_set.hdf5'),'/training_data_real');
received_imag = h5read(strcat('test_dataset_v3_',num2str(M),'_pilots_',num2str(SNR),'_data_set.hdf5'),'/training_data_imag');
received_signal = received_real + 1i*received_imag;
load(strcat('prec_comb_sig_',num2str(M),'_pilots_',num2str(SNR),'_data_set.mat'))

Phi = zeros(M*Lr,Nt*Nr); % Initialize measurement matrix

for ii = 1:M
    signal = signal_save(:,ii);
    Phi((ii-1)*Lr+(1:Lr),:)=kron(signal.'*Ftr_save(:,(ii-1)*Lt+(1:Lt)).',Wtr_save(:,(ii-1)*Lr+(1:Lr))');
end

Cwhiten = zeros(M*Lr,M*Lr);
for m = 1:M
    temporPart = Wtr_save(:,(m-1)*Lr+1:m*Lr);
    Cwhiten((m-1)*Lr+1:m*Lr,(m-1)*Lr+1:m*Lr) = temporPart'*temporPart;
end


RXanglesBar = linspace(0,pi,96*2);
TXanglesBar = linspace(0,pi,24*2);


RXanglesTilde = zeros(1,96*2);
TXanglesTilde = zeros(1,24*2);

Psi0 = dictionary_angles(RXanglesBar,TXanglesBar,RXanglesTilde,TXanglesTilde,Nr,Nt);

Psi0 = Psi0(:,vec(muuPowSel)>0);


for nc = 1:1000
    nc
    receivB = reshape(received_signal(nc,:,:),M*Lr,256);
    receivC = sqrt(256)*ifft(receivB,[],2);
    
    y = sqrtm(Cwhiten)\receivC;
    
    yPow = sum(abs(y).^2,1);
    [sortedd,indexxSort] = sort(yPow,'descend');
    
    K = max(5,length(find(sortedd>(max(sortedd)-median(sortedd))*0.05+median(sortedd))));
    tapIndex = indexxSort(1:K);
    y = y(:,tapIndex);
    
    gammaa = (gammaLow+gammaUpp)/2;
    
    
    Psi = (sqrtm(Cwhiten)\Phi)*Psi0;
    Psi2 = (Psi'*Psi);
    
    yy2 = Psi'*y;
    
    alpha = 0.1*ones(Gr*Gt,1);
    xOld = zeros(Gr*Gt,K);
    diffEM = 100;
    epsilonEM = 10^(-4);
    maxIterEM = 50;
    minIterME = 8;
    iterEM = 0;
    
    %Start iterations
    muu = zeros(Gr*Gt,K);
    chii  = repmat(1./alpha, [1,K]);
    
    while ((iterEM < maxIterEM) && (diffEM > epsilonEM)) || (iterEM<minIterME)
        
        etaa = alpha;
        for rr = find(vec(indicess(1:2:end,1:2:end)).'>0)
            etaa(rr) = etaa(rr) + betaa*sum(alpha(muuNeighbors(rr,muuNeighbors(rr,:)>0)));
        end
        
        temporInv = gammaa*Psi2+diag(etaa);
        chii = real(repmat(diag(temporInv\eye(Gr*Gt)),[1,K]));
        gammaDenom = sum((1-repmat(etaa,[1,K]).*chii)/gammaa,'all');
        muu=temporInv\(gammaa*yy2);
        
        
        gammaDenom = gammaDenom + sum(abs(y-Psi*muu).^2,'all');
        
        
        %M-Step: Hyper-parameter Learning
        
        temporVect = sum(abs(muu).^2 + chii,2);
        
        
        omegaa = temporVect;
        for rr = find(vec(indicess(1:2:end,1:2:end)).'>0)
            omegaa(rr) = omegaa(rr) + betaa*sum(temporVect(muuNeighbors(rr,muuNeighbors(rr,:)>0)));
        end
        
        
        alpha = K./(omegaa);
        
        
        %Step 2: Update for \gamma
        
        
        gammaa = (M*Lr*K)/gammaDenom;
        if gammaa < gammaLow
            gammaa = gammaLow;
        end
        if gammaa > gammaUpp
            gammaa = gammaUpp;
        end
        
        
        
        
        diffEM = sum(abs(xOld-muu).^2,'all')/sum(abs(muu).^2,'all');
        xOld = muu;
        
        iterEM  = iterEM+1;
        
    end
    muu2 = zeros(Gr*Gt,256);
    
    muu2(:,tapIndex)=muu;
    
    
    muu2(abs(muu2)<max(vec(abs(muu)))*0.08)=0;
    estimated_channels(nc,:,:,:) = reshape(fft(Psi0*muu2,256,2)/sqrt(256),Nr,Nt,256);
end

%%% save(strcat('channel_estimates_',num2str(M),'_pilots_',num2str(SNR),'_data_set.mat'),'estimated_channels','-v7.3');
