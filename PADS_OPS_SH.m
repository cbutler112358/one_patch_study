% One patch study of bi-sex harming homing drive, prototype of larger Patch
% Drive Simulator (PaDS). 
% Author: Cole Butler 
%
% 
%          ______   ________   ______   ______      
%         /_____/\ /_______/\ /_____/\ /_____/\     
%         \:::_ \ \\::: _  \ \\:::_ \ \\::::_\/_    
%          \:(_) \ \\::(_)  \ \\:\ \ \ \\:\/___/\   
%           \: ___\/ \:: __  \ \\:\ \ \ \\_::._\:\  
%            \ \ \    \:.\ \  \ \\:\/.:| | /____\:\ 
%             \_\/     \__\/\__\/ \____/_/ \_____\/ 
%                                           
% 
% 
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% A program simulating model dynamics on a local scale 
% (one patch, two adjacent patches, or two adjacent patches and a non-
% target patch). Simulation iterates daily and uses the life stages
% described on page J3;71. 
%
% INPUTS:
% -- NUM_GENS:
% -- NUM_GENS_RELEASE:
% -- driveParams:
% -- orgParams:
% -- dispParams:
% -- fitnessType: 
%
% OUTPUTS:
% -- popStruct: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [popStruct] = PADS_OPS_SH(NUM_GENS,NUM_GENS_RELEASE,...
    driveParams,orgParams,dispParams,fitnessType)

    filename = "SH_breeding_table.xlsx";
    raw = readcell(filename);

    % drive parameters
    s1 = driveParams(1);                 % fitness cost of homing drive;
    DOMINANCE = driveParams(2);         % dominance of fitness;
    CONV_EFFICIENCY = driveParams(3);   % conversion efficiency;
    RELEASE_RATIO = driveParams(4);     % how many transgenics are 
                                        % released?
    RELEASE_WINDOW = driveParams(5);    % if > 0, 3 multiple releases 
                                        % taking place in intervals of this
                                        % many days
    s2 = driveParams(6);                % fitness cost of Cas9/split 
                                        % construct
    if (RELEASE_WINDOW > 0)
        MULTI_RELEASE = true;
    else 
        MULTI_RELEASE = false;    
    end
                                        
    % mosquito parameters
    NU = orgParams(1);                  % daily density-independent 
                                        % survival probability for 
                                        % immatures (0.9351);
    ALPHA = orgParams(2);               % density-dependent parameter;
    BETA = orgParams(3);                % strength of density-dependence;
    LAMBDA = orgParams(4);              % no. of eggs per batch (86);
    MSR = orgParams(5);                 % male adult survival rate 
                                        % (0.6930);
    FSR = orgParams(6);                 % female adult survival rate 
                                        % (0.8249);
    OVIPROB = orgParams(7);             % probability of female ovipositing
                                        % (prob. of 0.2857 assumes avg 
                                        % gonotrophic cycle of 3.5 days);
    
    % migration and immigration probs
    MIG_PROB = dispParams(1);           % movement between local patches;
    IMMIG_PROB = dispParams(2);         % movement between global patches

    
    genoFreqs = raw(2:end,3:end); 
    % convert to string
    test = string(genoFreqs(:,:));
    % genotype conversion rate 
    genoFreqs = replace(test,'h',sprintf('%.2f',CONV_EFFICIENCY));
    % fitness cost that affects hatching probability
    % genoFreqs = replace(genoFreqs,'s',sprintf('%.2f',FITNESS_COST));
    
    % construct new array 
    zygoteFreqs = zeros(size(genoFreqs));
    matDim = size(genoFreqs);
    % zygote frequency matrix...
    for i = 1:matDim(1)
        for j = 1:matDim(2)
            zygoteFreqs(i,j) = eval(genoFreqs(i,j));
        end
    end

    % add a column for eggs that are not produced; zygoteFreqs follows male
    % first pairing specification
    deathProbVec = ones(length(zygoteFreqs),1) - sum(zygoteFreqs,2);
    zygoteFreqs = [zygoteFreqs, deathProbVec];

    %% define arrays for data storage
    
    NUM_GENOTYPES = 9;                          % no. of genotypes                        
    alleleFreqVec = zeros(1,NUM_GENS+1,3);      % allele frequency vector
    % [AACC, AACc, AAcc, 3
    % AaCC, AaCc, Aacc,  6
    % aaCC, aaCc, aacc]  9
    
    fitnessVec = [1, 1-DOMINANCE*s1, 1-s1, ...      % relative fitness vector
        (1-s2), (1-s2)*(1-DOMINANCE*s1), (1-s1)*(1-s2),...
        (1-s2)^2, (1-DOMINANCE*s1)*(1-s2)^2, (1-s1)*(1-s2)^2];      

    % matrices for storing male and female adults, etc. -- dimensions are
    % organized by generations, genotype, patch no.
    maleMat = zeros(NUM_GENS+1, NUM_GENOTYPES,3);
    femaleMat = zeros(NUM_GENS+1, NUM_GENOTYPES,3);
    maleEggMat = zeros(NUM_GENS+1, NUM_GENOTYPES,3);
    femaleEggMat = zeros(NUM_GENS+1, NUM_GENOTYPES,3);
    maleAquaticMat = zeros(NUM_GENS+1, NUM_GENOTYPES,3);
    femaleAquaticMat = zeros(NUM_GENS+1, NUM_GENOTYPES,3);
    % records mosquitoes during pre-adult life stages (10 days total)
    maleStageRecordMat = zeros(10, NUM_GENOTYPES,3);
    femaleStageRecordMat = zeros(10, NUM_GENOTYPES,3);
    totalStageRecordMat = maleStageRecordMat + femaleStageRecordMat; 
    
    % 0 generation, begin with INIT_POP male and female adults in each 
    % patch 
    INIT_POP = 100;
    maleMat(1,1,:) = INIT_POP;
    femaleMat(1,1,:) = INIT_POP;
    % initialize remaining matrices
    totalMat = maleMat + femaleMat;
    alleleFreqVec(1,1,:) = (totalMat(1,2,:) + 2*totalMat(1,3,:))./(2*sum(totalMat(1,:,:)));
    % matrix for females that have not oviposited
    femaleOvipositAllowed = zeros(NUM_GENOTYPES,NUM_GENOTYPES,3); 
    % [femaleMat(1,:,1); femaleMat(1,:,2); femaleMat(1,:,3)];
    % femaleTotalGeno = zeros(3,3,3);             % variable to keep track of 
                                                % females and the genotype
                                                % of the males they mate
                                                % with, regardless of
                                                % oviposition status
    RELEASE_COUNT = 0;
    %% run the sim
    % main for loop
    for i = 2:(NUM_GENS+1)
        % disp(i); 
        % do we release transgenics?
        if (i-1 == NUM_GENS_RELEASE)
            numRelease = RELEASE_RATIO*sum(maleMat(i-1,:,1));
            % release males
            maleMat(i-1,9,1) = ceil(numRelease);
            % no. of releases, all of the same size
            RELEASE_COUNT = 1;
        end
        % are we performing multiple releases?
        if (MULTI_RELEASE)
            if (RELEASE_COUNT < 3) && (i-1 == NUM_GENS_RELEASE + RELEASE_COUNT*RELEASE_WINDOW)  
                % release males
                maleMat(i-1,9,1) = ceil(numRelease);                
                RELEASE_COUNT = RELEASE_COUNT + 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%% AGING %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % density-dependent mortality for L3
        L3_pop = [sum(totalStageRecordMat(7,:,1)), sum(totalStageRecordMat(7,:,2)), ...
            sum(totalStageRecordMat(7,:,3))];
        survivalProbability = 1./(1+(ALPHA*L3_pop).^BETA);
        maleStageRecordMat(7,:,1) = binornd(maleStageRecordMat(7,:,1),survivalProbability(1));
        femaleStageRecordMat(7,:,1) = binornd(femaleStageRecordMat(7,:,1),survivalProbability(1));
        maleStageRecordMat(7,:,2) = binornd(maleStageRecordMat(7,:,2),survivalProbability(2));
        femaleStageRecordMat(7,:,2) = binornd(femaleStageRecordMat(7,:,2),survivalProbability(2));
        maleStageRecordMat(7,:,3) = binornd(maleStageRecordMat(7,:,3),survivalProbability(3));
        femaleStageRecordMat(7,:,3) = binornd(femaleStageRecordMat(7,:,3),survivalProbability(3));        
        
        % before transitions, apply fitness cost to relevant life stage
        if (fitnessType == "EA")
            % early-acting fitness cost imposed b/w egg and L1, i.e. row 3
            maleStageRecordMat(3,:,1) = binornd(maleStageRecordMat(3,:,1), fitnessVec);
            femaleStageRecordMat(3,:,1) = binornd(femaleStageRecordMat(3,:,1), fitnessVec);
            maleStageRecordMat(3,:,2) = binornd(maleStageRecordMat(3,:,2), fitnessVec);
            femaleStageRecordMat(3,:,2) = binornd(femaleStageRecordMat(3,:,2), fitnessVec);
            maleStageRecordMat(3,:,3) = binornd(maleStageRecordMat(3,:,3), fitnessVec);
            femaleStageRecordMat(3,:,3) = binornd(femaleStageRecordMat(3,:,3), fitnessVec);
        end

        if (fitnessType == "LA")
            % late-acting fitness cost imposed b/w pupa and adult, i.e. row
            % 10
            maleStageRecordMat(10,:,1) = binornd(maleStageRecordMat(10,:,1), fitnessVec);
            femaleStageRecordMat(10,:,1) = binornd(femaleStageRecordMat(10,:,1), fitnessVec);
            maleStageRecordMat(10,:,2) = binornd(maleStageRecordMat(10,:,2), fitnessVec);
            femaleStageRecordMat(10,:,2) = binornd(femaleStageRecordMat(10,:,2), fitnessVec);
            maleStageRecordMat(10,:,3) = binornd(maleStageRecordMat(10,:,3), fitnessVec);
            femaleStageRecordMat(10,:,3) = binornd(femaleStageRecordMat(10,:,3), fitnessVec);
        end        
        
        % advance a life stage for each immature...
        emergingMaleAdults = reshape(maleStageRecordMat(10,:,:),NUM_GENOTYPES,3);
        emergingFemaleAdults = reshape(femaleStageRecordMat(10,:,:),NUM_GENOTYPES,3);
        % add newly emerged adults to each patch
        maleMat(i,:,:) = reshape(emergingMaleAdults,1,NUM_GENOTYPES,3);
        femaleMat(i,:,:) = reshape(emergingFemaleAdults,1,NUM_GENOTYPES,3);
        % only newly emerged females can mate
        if (i == 2)
            % no eggs upon initialization
            femaleBreedAllowed = reshape(femaleMat(i-1,:,:),NUM_GENOTYPES,3)';
        else
            femaleBreedAllowed = reshape(femaleMat(i,:,:),NUM_GENOTYPES,3)';
        end
        % females that have not yet oviposited
        %%% femaleOvipositAllowed = femaleOvipositAllowed + femaleBreedAllowed;
        % shift all the stages one down
        maleStageRecordMat(2:end,:,:) = maleStageRecordMat(1:9,:,:);
        maleStageRecordMat(1,:,:) = 0;
        femaleStageRecordMat(2:end,:,:) = femaleStageRecordMat(1:9,:,:);
        femaleStageRecordMat(1,:,:) = 0;    
        
        % updated counts are (i) adults that emerged and (ii) adults
        % surviving previous day 
        maleMat(i,:,:) = maleMat(i,:,:) + maleMat(i-1,:,:);
        femaleMat(i,:,:) = femaleMat(i,:,:) + femaleMat(i-1,:,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%% MATING %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % fprintf("Simulating generation %.f of %.f\n",i-1,NUM_GENS);

        % ----------------------
        %%%%%%% PATCH 1 %%%%%%%%
        % ----------------------

        % probability of mating with male of a certain genotype
        maleProbVec_p1 = maleMat(i,:,1)/sum(maleMat(i,:,1));
        % matings of females with males; female genotype per column, male
        % genotype per row
        malePairings = mnrnd(femaleBreedAllowed(1,:)', maleProbVec_p1)';
        % if there are no males present, there are zero pairings (MATLAB
        % returns NaNs)
        malePairings(isnan(malePairings)) = 0;
        % femaleTotalGeno(:,:,1) = femaleTotalGeno(:,:,1) + malePairings; 
        femaleOvipositAllowed(:,:,1) = femaleOvipositAllowed(:,:,1) + malePairings;
        % females and males disperse before oviposition
        female_p1_to_p2 = binornd(femaleOvipositAllowed(:,:,1),MIG_PROB);
        femaleOvipositAllowed(:,:,1) = femaleOvipositAllowed(:,:,1) - female_p1_to_p2;
        male_p1_to_p2 = binornd(maleMat(i,:,1),MIG_PROB); 
        maleMat(i,:,1) = maleMat(i,:,1) - male_p1_to_p2; 

        % ----------------------
        %%%%%%% PATCH 2 %%%%%%%%
        % ----------------------            

        maleProbVec_p2 = maleMat(i,:,2)/sum(maleMat(i,:,2));
        malePairings = mnrnd(femaleBreedAllowed(2,:)', maleProbVec_p2)';
        malePairings(isnan(malePairings)) = 0;
        % femaleTotalGeno(:,:,2) = femaleTotalGeno(:,:,2) + malePairings;
        femaleOvipositAllowed(:,:,2) = femaleOvipositAllowed(:,:,2) + malePairings;
        female_p2_to_p1 = binornd(femaleOvipositAllowed(:,:,2),MIG_PROB);
        femaleOvipositAllowed(:,:,2) = femaleOvipositAllowed(:,:,2) - female_p2_to_p1;
        maleDispersal = binornd(maleMat(i,:,2),MIG_PROB + IMMIG_PROB);
        maleMat(i,:,2) = maleMat(i,:,2) - maleDispersal;
        if MIG_PROB + IMMIG_PROB > 0
            male_p2_to_p1 = binornd(maleDispersal, MIG_PROB/(MIG_PROB + IMMIG_PROB));
        else
            % both probabilities are zero
            male_p2_to_p1 = binornd(maleDispersal, MIG_PROB);
        end
        % dispersal (immigration) to patch 3 as well...
        female_p2_to_p3 = binornd(femaleOvipositAllowed(:,:,2),IMMIG_PROB);
        femaleOvipositAllowed(:,:,2) = femaleOvipositAllowed(:,:,2) - female_p2_to_p3;
        male_p2_to_p3 = maleDispersal - male_p2_to_p1;       

        % ----------------------
        %%%%%%% PATCH 3 %%%%%%%%
        % ---------------------- 

        maleProbVec_p3 = maleMat(i,:,3)/sum(maleMat(i,:,3));
        malePairings = mnrnd(femaleBreedAllowed(3,:)', maleProbVec_p3)';
        malePairings(isnan(malePairings)) = 0;
        % femaleTotalGeno(:,:,3) = femaleTotalGeno(:,:,3) + malePairings;
        femaleOvipositAllowed(:,:,3) = femaleOvipositAllowed(:,:,3) + malePairings;
        female_p3_to_p2 = binornd(femaleOvipositAllowed(:,:,3),IMMIG_PROB);
        femaleOvipositAllowed(:,:,3) = femaleOvipositAllowed(:,:,3) - female_p3_to_p2;
        male_p3_to_p2 = binornd(maleMat(i,:,3), IMMIG_PROB);
        maleMat(i,:,3) = maleMat(i,:,3) - male_p3_to_p2;  

        %%% update adult composition in all patches
        femaleOvipositAllowed(:,:,1) = femaleOvipositAllowed(:,:,1) + female_p2_to_p1;
        femaleOvipositAllowed(:,:,2) = femaleOvipositAllowed(:,:,2) + female_p1_to_p2 + female_p3_to_p2;
        femaleOvipositAllowed(:,:,3) = femaleOvipositAllowed(:,:,3) + female_p2_to_p3;
        
        % femaleTotalGeno(:,:,1) = femaleTotalGeno(:,:,1) - female_p1_to_p2 + female_p2_to_p1;
        % femaleTotalGeno(:,:,2) = femaleTotalGeno(:,:,2) + female_p1_to_p2 + female_p3_to_p2 - female_p2_to_p1 - female_p2_to_p3;
        % femaleTotalGeno(:,:,3) = femaleTotalGeno(:,:,3) + female_p2_to_p3 - female_p3_to_p2;
        
        femaleMat(i,:,1) = femaleMat(i,:,1) - sum(female_p1_to_p2) + sum(female_p2_to_p1);
        femaleMat(i,:,2) = femaleMat(i,:,2) + sum(female_p1_to_p2) + sum(female_p3_to_p2) - sum(female_p2_to_p1) - sum(female_p2_to_p3);
        femaleMat(i,:,3) = femaleMat(i,:,3) + sum(female_p2_to_p3) - sum(female_p3_to_p2);
        
        maleMat(i,:,1) = maleMat(i,:,1) + male_p2_to_p1; %  - male_p1_to_p2;
        maleMat(i,:,2) = maleMat(i,:,2) + male_p1_to_p2 + male_p3_to_p2; %  - male_p2_to_p3 - male_p2_to_p1;
        maleMat(i,:,3) = maleMat(i,:,3) + male_p2_to_p3; %  - male_p3_to_p2;
           

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% OVIPOSITION %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

        % no. of eggs produced from pairings (oviposition occurs
        % probabilistically)
        ovipositFemaleMat_p1 = binornd(femaleOvipositAllowed(:,:,1),OVIPROB);
        ovipositFemaleMat_p2 = binornd(femaleOvipositAllowed(:,:,2),OVIPROB);
        ovipositFemaleMat_p3 = binornd(femaleOvipositAllowed(:,:,3),OVIPROB);
        numOffspringMat_p1 = poissrnd(ovipositFemaleMat_p1*LAMBDA);
        numOffspringMat_p2 = poissrnd(ovipositFemaleMat_p2*LAMBDA);
        numOffspringMat_p3 = poissrnd(ovipositFemaleMat_p3*LAMBDA);

        % number of eggs by zygote stored in numOffspringVec
        % (organized by male genotype first)
        pairingVec = reshape(numOffspringMat_p1',1,[]);
        offspringGenotypeVec = sum(mnrnd(pairingVec',zygoteFreqs));
        % remove eggs that are not produced
        offspringGenotypeVec = offspringGenotypeVec(1:3); 
        % determine no. of female and male eggs produced from mating
        maleEggs_p1 = binornd(offspringGenotypeVec,0.5);
        femaleEggs_p1 = offspringGenotypeVec - maleEggs_p1; 

        % do the same for each patch...
        pairingVec = reshape(numOffspringMat_p2',1,[]);
        offspringGenotypeVec = sum(mnrnd(pairingVec',zygoteFreqs));
        offspringGenotypeVec = offspringGenotypeVec(1:3); 
        maleEggs_p2 = binornd(offspringGenotypeVec,0.5);
        femaleEggs_p2 = offspringGenotypeVec - maleEggs_p2; 

        pairingVec = reshape(numOffspringMat_p3',1,[]);
        offspringGenotypeVec = sum(mnrnd(pairingVec',zygoteFreqs));
        offspringGenotypeVec = offspringGenotypeVec(1:3); 
        maleEggs_p3 = binornd(offspringGenotypeVec,0.5);
        femaleEggs_p3 = offspringGenotypeVec - maleEggs_p3; 

        % update stage recorders
        maleStageRecordMat(1,:,:) = reshape([maleEggs_p1; maleEggs_p2; maleEggs_p3]',1,NUM_GENOTYPES,3);
        femaleStageRecordMat(1,:,:) = reshape([femaleEggs_p1; femaleEggs_p2; femaleEggs_p3]',1,NUM_GENOTYPES,3);
        
        % females assumed to oviposit all at once, remove these from the
        % matrix containing those who have not oviposited
        femaleOvipositAllowed(:,:,1) = femaleOvipositAllowed(:,:,1) - ovipositFemaleMat_p1;
        femaleOvipositAllowed(:,:,2) = femaleOvipositAllowed(:,:,2) - ovipositFemaleMat_p2;
        femaleOvipositAllowed(:,:,3) = femaleOvipositAllowed(:,:,3) - ovipositFemaleMat_p3;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%% DEATH %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

        % immature density-independent death
        maleStageRecordMat = binornd(maleStageRecordMat, NU); 
        femaleStageRecordMat = binornd(femaleStageRecordMat, NU);
        
        % adult death        
        % for females, we need to distinguish between those that have
        % oviposited and those that have not--first, find the females that
        % cannot oviposit and kill them
        survivingNonOvipositFemales = binornd(femaleMat(i,:,:) - sum(femaleOvipositAllowed),FSR);
        % death of female adults who can oviposit
        femaleOvipositAllowed = binornd(femaleOvipositAllowed, FSR); 
        % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %%% this should ROUGHLY be binornd(femaleMat(i,:,:),FSR) %%%
        % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        femaleMat(i,:,:) = survivingNonOvipositFemales + sum(femaleOvipositAllowed); 
        % death of male adults
        maleMat(i,:,:) = binornd(maleMat(i,:,:),MSR); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% UPDATE ARRAYS %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        %%% update data matrices
        maleEggMat(i,:,:) = sum(maleStageRecordMat(1:3,:,:));
        femaleEggMat(i,:,:) = sum(femaleStageRecordMat(1:3,:,:));
        maleAquaticMat(i,:,:) = sum(maleStageRecordMat);
        femaleAquaticMat(i,:,:) = sum(femaleStageRecordMat);
        totalMat(i,:,:) = femaleMat(i,:,:) + maleMat(i,:,:); 
        alleleFreqVec(1,i,:) = (totalMat(i,2,:) + 2*totalMat(i,3,:))./(2*sum(totalMat(i,:,:)));
        
        % update total counts of immature stages
        totalStageRecordMat = maleStageRecordMat + femaleStageRecordMat; 
    end
    
    popStruct = struct();
    popStruct.femaleMat = femaleMat;
    popStruct.maleMat = maleMat;
    popStruct.femaleAquaticMat = femaleAquaticMat;
    popStruct.maleAquaticMat = femaleAquaticMat;
    
    return
end





