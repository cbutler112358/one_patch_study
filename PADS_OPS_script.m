

close all 

NUM_GENS = 365;
NUM_GENS_RELEASE = 100;
driveParams = [0.66, 0.83, 0.95, 6, 7];
orgParams = [0.9351, 0.001, 2, 86, 0.6930, 0.8249, 0.2857];
dispParams = [0.1305, 0.01];
fitnessType = 'LA';

tmp = PADS_OPS(NUM_GENS,NUM_GENS_RELEASE,driveParams,orgParams,dispParams,...
    fitnessType);


% plot males and females from each cell
subplot(1,3,1)
plot(tmp.femaleMat(:,1,1),'-k','linewidth',2);
hold on
plot(tmp.femaleMat(:,2,1),'-b','linewidth',2);
plot(tmp.femaleMat(:,3,1),'-r','linewidth',2);

subplot(1,3,2)
plot(tmp.femaleMat(:,1,2),'-k','linewidth',2);
hold on
plot(tmp.femaleMat(:,2,2),'-b','linewidth',2);
plot(tmp.femaleMat(:,3,2),'-r','linewidth',2);

subplot(1,3,3)
plot(tmp.femaleMat(:,1,3),'-k','linewidth',2);
hold on
plot(tmp.femaleMat(:,2,3),'-b','linewidth',2);
plot(tmp.femaleMat(:,3,3),'-r','linewidth',2);

figure 
subplot(1,3,1)
plot(tmp.maleMat(:,1,1),'-k','linewidth',2);
hold on
plot(tmp.maleMat(:,2,1),'-b','linewidth',2);
plot(tmp.maleMat(:,3,1),'-r','linewidth',2);

subplot(1,3,2)
plot(tmp.maleMat(:,1,2),'-k','linewidth',2);
hold on
plot(tmp.maleMat(:,2,2),'-b','linewidth',2);
plot(tmp.maleMat(:,3,2),'-r','linewidth',2);

subplot(1,3,3)
plot(tmp.maleMat(:,1,3),'-k','linewidth',2);
hold on
plot(tmp.maleMat(:,2,3),'-b','linewidth',2);
plot(tmp.maleMat(:,3,3),'-r','linewidth',2);
