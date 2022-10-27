% Batch simulation studying the invasion potential of a simple homing
% drive. In this script, we study how a simple homing drive for varying
% values of h and s can possibly invade Patch 2 (and possibly Patch 3) if
% Patch 1 has been invaded.

migVec = 0:0.02:0.5;
immigVec = 0:0.0005:0.01;

% Matrices to store everything in!
alleleFreqMat_p1 = zeros(length(hVec), length(sVec));
alleleFreqMat_p2 = zeros(length(hVec), length(sVec));
alleleFreqMat_p3 = zeros(length(hVec), length(sVec));
meanPopMat_p1 = zeros(length(hVec), length(sVec));
meanPopMat_p2 = zeros(length(hVec), length(sVec));
meanPopMat_p3 = zeros(length(hVec), length(sVec));

% The fixed population is allowed to reach equilibrium for 200 days.
% After this, immigration is "turned on" between Patch 2 and Patch 3 and
% the simulation allowed to run for 1 year. Invasion of Patch 3 occurs if 
% transgene frequency exceeds 5%. 

NUM_GENS = 565+365;
NUM_GENS_RELEASE = 200;
driveParams = [0.74, 0.68, 0.95, 0, 16, 0];
orgParams = [0.9351, 0.001, 2, 86, 0.6930, 0.8249, 0.2857];
fitnessType = 'LA';
releaseInd = 2;
homoInd = [1,1];

for i = 1:length(migVec)
    for j = 1:length(immigVec)
        disp([i,j]); 
        dispParams = [migVec(i), immigVec(j)];


        tmp = PADS_OPS(NUM_GENS,NUM_GENS_RELEASE,driveParams,orgParams,dispParams,...
            fitnessType, releaseInd, homoInd);        
        
        % calculate transgene frequencies in each patch
        totalMat = tmp.maleMat + tmp.femaleMat;
        
        % Store (i) allele frequency and (ii) population per patch. The
        % population will be calculated as a an average over the final 100 
        % days (after immigration is "turned on"), while allele frequency 
        % will be considered as a snapshot.
        
        alleleFreq_p1 = (2*totalMat(end,1,1) + totalMat(end,2,1))./(2*sum(totalMat(end,:,1),2));
        alleleFreq_p2 = (2*totalMat(end,1,2) + totalMat(end,2,2))./(2*sum(totalMat(end,:,2),2));
        alleleFreq_p3 = (2*totalMat(end,1,3) + totalMat(end,2,3))./(2*sum(totalMat(end,:,3),2));
        
        meanPop_p1 = mean(sum(totalMat((end-99):end,:,1),2));
        meanPop_p2 = mean(sum(totalMat((end-99):end,:,2),2));
        meanPop_p3 = mean(sum(totalMat((end-99):end,:,3),2));
        
        % Average WT pop. (Patch 3)
        avgEquilWT_pop = mean(sum(totalMat(150:199,:,3),2));
        
        % Store data points.
        alleleFreqMat_p1(i,j) = alleleFreq_p1;
        alleleFreqMat_p2(i,j) = alleleFreq_p2;
        alleleFreqMat_p3(i,j) = alleleFreq_p3;
        meanPopMat_p1(i,j) = meanPop_p1/avgEquilWT_pop;
        meanPopMat_p2(i,j) = meanPop_p2/avgEquilWT_pop;
        meanPopMat_p3(i,j) = meanPop_p3/avgEquilWT_pop;
    end
end