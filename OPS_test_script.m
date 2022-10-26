% Just a general script to test the performance of various functions
% related to the one patch study files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% some code simulating the Drosophila system developed by Max. Uses a
% sex-linked tethered homing drive. This was just testing whether automatic
% generation of the breeding table worked as expected.
rho = 1;
MALE_CONV_RATE = 0.4136;
FEMALE_CONV_RATE = 0.2262;
% MALE_CONV_RATE = 0.1382;
% FEMALE_CONV_RATE = NaN;
multiRelease = true;
fitnessCostVec = [0.05, 0.05];
RELATIVE_FECUNDITY = 0.727;
graphBool = true;

dataMat = cage_trial_split_rawBT(multiRelease,rho,...
    MALE_CONV_RATE,FEMALE_CONV_RATE,fitnessCostVec,RELATIVE_FECUNDITY,graphBool);

test = genBT_simple_tethered_homing(MALE_CONV_RATE,FEMALE_CONV_RATE, ...
    RELATIVE_FECUNDITY,1);


dataMat2 = cage_trial_split(multiRelease,rho,...
    MALE_CONV_RATE,FEMALE_CONV_RATE,fitnessCostVec,RELATIVE_FECUNDITY,graphBool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% More code testing automatic breeding table generation.
test = genBT_TLU_tethered_homing(MALE_CONVERSION_PROB,FEMALE_CONVERSION_PROB);

test = genBT_tethered_homing(MALE_CONVERSION_PROB,FEMALE_CONVERSION_PROB);

test = genBT_simple_homing(CONV_EFFICIENCY,CONV_EFFICIENCY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Seeing how the frequency of a transgene and it's "anchor" evolve with 
% time in the one putch study for a tethered homing drive.
NUM_GENS = 1250;
NUM_GENS_RELEASE = 100;
driveParams = [0.6, 0.75, 1.0, 0, 8, 14];
orgParams = [0.9351, 0.001, 2, 86, 0.6930, 0.8249, 0.2857];
dispParams = [0,0]; % [0.1305, 0.001];
fitnessType = 'LA';
releaseInd = 3;

tmp = PADS_OPS(NUM_GENS,NUM_GENS_RELEASE,driveParams,orgParams,dispParams,...
    fitnessType, releaseInd);

% allele frequency vectors
totalMat = tmp.maleMat + tmp.femaleMat;
cFreqVec = (2*sum(totalMat(:,1:3,1),2) + sum(totalMat(:,4:6,1),2))./(2*sum(totalMat(:,:,1),2));
gRNAFreqVec = (2*sum(totalMat(:,[1,4,7],1),2) + sum(totalMat(:,[2,5,8],1),2))./(2*sum(totalMat(:,:,1),2));

% plot(sum(tmp.maleMat(:,1:8,1),2) + sum(tmp.femaleMat(:,1:8,1),2),'-r')
% hold on
% plot(tmp.maleMat(:,end,1) + tmp.femaleMat(:,end,1),'-b')
subplot(1,2,1)
plot(cFreqVec,'-r')
ylim([0,1])
subplot(1,2,2)
plot(gRNAFreqVec,'-r')
ylim([0,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some code testing a new implementation in the OPS code. Specifically, I 
% realized that simulating a full three patch system was not necessary to
% study the potential of the drive spreading from one patch to another.
% That being the case, the new "homoInd" variable initializes patches 1 and
% 2 with a pop fixed for that index. Plots are an added bonus...
NUM_GENS = 750;
NUM_GENS_RELEASE = 150;
driveParams = [0.8, 0.65, 1.0, 0, 16, 0];
orgParams = [0.9351, 0.001, 2, 86, 0.6930, 0.8249, 0.2857];
dispParams = [0.1305, 0.001];
fitnessType = 'LA';
releaseInd = 2;
homoInd = -1;

tmp = PADS_OPS(NUM_GENS,NUM_GENS_RELEASE,driveParams,orgParams,dispParams,...
    fitnessType, releaseInd, homoInd);

% allele frequency vectors
totalMat = tmp.maleMat + tmp.femaleMat;
% cFreqVec = (2*sum(totalMat(:,1:3,1),2) + sum(totalMat(:,4:6,1),2))./(2*sum(totalMat(:,:,1),2));
alleleFreqVec = (2*sum(totalMat(:,1,1),2) + sum(totalMat(:,2,1),2))./(2*sum(totalMat(:,:,1),2));
subplot(1,3,1)
plot(alleleFreqVec,'-r');
ylim([0,1]);
alleleFreqVec = (2*sum(totalMat(:,1,2),2) + sum(totalMat(:,2,2),2))./(2*sum(totalMat(:,:,2),2));
subplot(1,3,2)
plot(alleleFreqVec,'-r');
ylim([0,1]);
alleleFreqVec = (2*sum(totalMat(:,1,3),2) + sum(totalMat(:,2,3),2))./(2*sum(totalMat(:,:,3),2));
subplot(1,3,3)
plot(alleleFreqVec,'-r');
ylim([0,1]);

figure 
% plot populations
subplot(1,3,1)
plot(sum(tmp.maleMat(:,:,1),2),'-b','linewidth',2);
hold on
plot(sum(tmp.femaleMat(:,:,1),2),'-m','linewidth',2);
subplot(1,3,2)
plot(sum(tmp.maleMat(:,:,2),2),'-b','linewidth',2);
hold on
plot(sum(tmp.femaleMat(:,:,2),2),'-m','linewidth',2);
subplot(1,3,3)
plot(sum(tmp.maleMat(:,:,3),2),'-b','linewidth',2);
hold on
plot(sum(tmp.femaleMat(:,:,3),2),'-m','linewidth',2);








