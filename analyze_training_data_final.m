muuPow = zeros(96,24);

for sett = 1:10
    
    load(strcat('muu_training_',num2str(sett),'.mat'));
    muuVec = reshape(muuTot,96,24,89*1000);
    muuPow = muuPow + sum(abs(muuVec).^2,3);
    
end

InterpFacRx = 2;
InterpFacTx = 2;

muuPowInt = zeros(96*InterpFacRx,24*InterpFacTx);
muuPowInt(1:InterpFacRx:end,1:InterpFacTx:end) = muuPow;

for rr = 1:96
    if rr<96
        for ii = 1:InterpFacRx-1
            muuPowInt((rr-1)*InterpFacRx+ii+1,:) = muuPowInt((rr-1)*InterpFacRx+1,:)+(muuPowInt(rr*InterpFacRx+1,:)- muuPowInt((rr-1)*InterpFacRx+1,:))/InterpFacRx*ii;
        end
    else
        for ii = 1:InterpFacRx-1
            muuPowInt((rr-1)*InterpFacRx+ii+1,:) = muuPowInt((rr-1)*InterpFacRx+1,:)+(0- muuPowInt((rr-1)*InterpFacRx+1,:))/InterpFacRx*ii;
        end
        
    end
    
    
end

for rr = 1:24
    if rr<24
        for ii = 1:InterpFacTx-1
            muuPowInt(:,(rr-1)*InterpFacTx+ii+1) = muuPowInt(:,(rr-1)*InterpFacTx+1)+(muuPowInt(:,rr*InterpFacTx+1)- muuPowInt(:,(rr-1)*InterpFacTx+1))/InterpFacTx*ii;
        end
    else
        for ii = 1:InterpFacTx-1
            muuPowInt(:,(rr-1)*InterpFacTx+ii+1) = muuPowInt(:,(rr-1)*InterpFacTx+1)+(0- muuPowInt(:,(rr-1)*InterpFacTx+1))/InterpFacTx*ii;
        end
        
    end
    
    
end


muuPowSel = zeros(size(muuPowInt));
muuPowSel(1:InterpFacRx:end,1:InterpFacTx:end) = 1;
[sortedd, ~]= sort(vec(muuPowInt),'descend');
iterr = 0;
while sum(muuPowSel,'all')<96*32
    iterr = iterr+1;
    threshold = sortedd(1+iterr);
    muuPowSel = muuPowInt>=threshold;
    muuPowSel(1:InterpFacRx:end,1:InterpFacTx:end) = 1;
end

maxDistVec = 7:-1:2;
intervalss= linspace(log(min(vec(muuPowInt))),log(max(vec(muuPowInt))),6);

stopp = 0;
iterr= 0;
while stopp<1
    iterr = iterr +1;
    maxx = 0;
    minn = inf;
    
    for rr = 1:96*InterpFacRx
        for tt = 1:24*InterpFacTx
            xIndices = max(1,rr-1):min(96*InterpFacRx,rr+1);
            yIndices = max(1,tt-1):min(24*InterpFacTx,tt+1);
            metricc = sum(muuPowInt(xIndices,yIndices),'all')/(length(xIndices)*length(yIndices));
            maxDist = maxDistVec(sum(intervalss<=log(metricc)));
            if (metricc>maxx) && (muuPowSel(rr,tt)<1)
                
                maxx = metricc;
                maxXindex = rr;
                maxYindex = tt;
            end
            if (metricc<minn) && (muuPowSel(rr,tt)>0)
                 xIndices2 = max(1,rr-maxDist):min(96*InterpFacRx,rr+maxDist);
                 yIndices2 = max(1,tt-maxDist):min(24*InterpFacTx,tt+maxDist);
                 conditX = 1;
                 conditY = 1;
                 
                 muuPowSelCandiX = muuPowSel(:,tt);
                 muuPowSelCandiY = muuPowSel(rr,:);
                 muuPowSelCandiX(rr) = 0;
                 muuPowSelCandiY(tt) = 0;
                 for rr2 = 1:length(xIndices2)-maxDist
                     
                    if sum(muuPowSelCandiX(xIndices2(rr2:rr2+maxDist)))<2
                        conditX = 0;
                    end
                 end
                 for tt2 = 1:length(yIndices2)-maxDist
                    if sum(muuPowSelCandiY(yIndices2(tt2:tt2+maxDist)))<2
                        conditY = 0;
                    end
                 end
                 if (conditX>0)&&(conditY>0)     
                     minn = metricc;
                     minXindex = rr;
                     minYindex = tt;
                 end
            end
        end
    end
    
    if minn > maxx
        break
    else
        muuPowSel(minXindex,minYindex) = 0;
        muuPowSel(maxXindex,maxYindex) = 1;

    end
    
end
               
            
            







indicess = zeros(size(muuPowSel));
muuNeighbors = zeros(96*32,4);
distancee = 3;
iterr = 0;
for tt = 1:24*InterpFacTx
    for rr = 1:96*InterpFacRx
        if muuPowSel(rr,tt)>0
            iterr = iterr+1;
            indicess(rr,tt) = iterr;
        end
    end
end

for tt = 1:24*InterpFacTx
    for rr = 1:96*InterpFacRx
        if muuPowSel(rr,tt)>0
            whichIndexx = indicess(rr,tt);
            startt = 0;
            for ii = 1:distancee-1
                if rr-ii>0
                    if muuPowSel(rr-ii,tt)>0
                        
                        startt = startt +1;
                
                        muuNeighbors(whichIndexx,startt) = indicess(rr-ii,tt);
                        break
                    end
                end
            end
            
            for ii = 1:distancee-1
                if rr+ii<=96*2
                    if muuPowSel(rr+ii,tt)>0

                        startt = startt +1;
                
                        muuNeighbors(whichIndexx,startt) = indicess(rr+ii,tt);
                        break
                    end
                end
            end
            for ii = 1:distancee-1
                if tt-ii>0
                    if muuPowSel(rr,tt-ii)>0

                        startt = startt +1;
                
                        muuNeighbors(whichIndexx,startt) = indicess(rr,tt-ii);
                        break
                    end
                end
            end
            
            for ii = 1:distancee-1
                if tt+ii<=24*2
                    if muuPowSel(rr,tt+ii)>0

                        startt = startt +1;
                
                        muuNeighbors(whichIndexx,startt) = indicess(rr,tt+ii);
                        break
                    end
                end
            end
        end
    end
end

save 'learned_parameters.mat' muuPowSel muuNeighbors indicess