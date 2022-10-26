% trial batch simulation for PADS_OPS (simple homing)


hVec = 0:0.1:1;
sVec = 0:0.1:1;
% no. of replicates for each set of conditions
REPLICATES = 5;

% matrices storing final pop. and allele freq.
finalPop = zeros(length(hVec), length(sVec), 3);
alleleFreq = zeros(length(hVec), length(sVec), 3);

NUM_GENS = 365*5;
NUM_GENS_RELEASE = 365;
orgParams = [0.9351, 0.001, 2, 86, 0.6930, 0.8249, 0.2857];
dispParams = [0.1305, 0.01];
fitnessType = 'LA';

for i = 1:length(hVec)
    for j = 1:length(sVec)
        disp([i,j]); 
        finalPopVec = zeros(1,REPLICATES,3);
        alleleFreqVec = zeros(1,REPLICATES,3);
        for k = 1:REPLICATES
            driveParams = [sVec(i), hVec(j), 0.95, 2, 14];
            
            tmp = PADS_OPS(NUM_GENS,NUM_GENS_RELEASE,driveParams,orgParams,dispParams,...
                fitnessType);
            
            totalMat = tmp.maleMat + tmp.femaleMat;
            totalMat = sum(totalMat(end,:,:),1);
            
            % first patch
            finalPopVec(1,k,1) = sum(totalMat(:,:,1));
            alleleFreqVec(1,k,1) = (totalMat(:,2,1) + 2*totalMat(:,3,1))/(2*finalPopVec(1,k,1));
            % second patch
            finalPopVec(1,k,2) = sum(totalMat(:,:,2));
            alleleFreqVec(1,k,2) = (totalMat(:,2,2) + 2*totalMat(:,3,2))/(2*finalPopVec(1,k,2));
            % third patch
            finalPopVec(1,k,3) = sum(totalMat(:,:,3));
            alleleFreqVec(1,k,3) = (totalMat(:,2,3) + 2*totalMat(:,3,3))/(2*finalPopVec(1,k,3));            
        end
        
        % first patch
        finalPop(i,j,1) = mean(finalPopVec(1,:,1));
        alleleFreq(i,j,1) = mean(alleleFreqVec(1,:,1));
        % second patch
        finalPop(i,j,2) = mean(finalPopVec(1,:,2));
        alleleFreq(i,j,2) = mean(alleleFreqVec(1,:,2));
        % third patch
        finalPop(i,j,3) = mean(finalPopVec(1,:,3));
        alleleFreq(i,j,3) = mean(alleleFreqVec(1,:,3));
    end
end