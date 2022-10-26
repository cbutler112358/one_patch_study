% Function algorithmically producing breeding tables for a simple homing 
% drive.
% Author: Cole Butler 
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The following function produces a breeding table given conversion
% efficiencies for either sex. Constructs (transgenes) are indicated by
% capital letters while wild-type are denoted by lower-case letters.
%
% Note that sex determination and fitness costs are incorporated externally
% (i.e. not within this function).
%
% INPUTS:
% -- MALE_CONVERSION_PROB: Conversion efficiency of construct in males.
% -- FEMALE_CONVERSION_PROB: Conversion efficiency of construct in females.
%
% OUTPUTS:
% -- breedingTable: A structure array containing the breeding table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [breedingTables] = genBT_simple_homing(MALE_CONVERSION_PROB, ...
    FEMALE_CONVERSION_PROB)

%% produce breeding table for male zygotes
tmpArray = readcell('SH_breeding_table.xlsx');

% pull out the part of the breeding table we will update
genoArray = tmpArray(2:end,3:end);
arrayDim = size(genoArray);
numPairings = arrayDim(1);
numGenotypes = arrayDim(2);
% include a column for nonviable progeny
genoArray = zeros(numPairings, numGenotypes); 

% dominant sterile--hemizygous males fertile, homozygotes and 
% hemizygous females sterile

% for each possible pairing, determine zygote frequencies...
for i = 1:numPairings
    % initialize haplotypes @ 0 (for male progeny, the father never
    % contributes a X chromosome)
    fatherHaplo = struct();    
    fatherHaplo.A = 0;
    fatherHaplo.a = 0;
    
    motherHaplo = fatherHaplo;    
    
    dadGenotype = string(tmpArray(1+i,1));
    momGenotype = string(tmpArray(1+i,2));

    % progeny are produced, check for homing in either parent
    homingDadBool = (dadGenotype == "Aa");
    homingMomBool = (momGenotype == "Aa");

    if (homingDadBool)
        % with MALE_CONVERSION_PROB, then father's genotype is CC
        fatherHaplo.A = (1-MALE_CONVERSION_PROB)*0.5 + MALE_CONVERSION_PROB;
        fatherHaplo.a = (1-MALE_CONVERSION_PROB)*0.5;
    else 
        % homing is not possible
        fatherHaplo.A = fatherHaplo.A + count(dadGenotype,"A")/2;
        fatherHaplo.a = fatherHaplo.a + count(dadGenotype,"a")/2;
    end
    
    if (homingMomBool)
        % with FEMALE_CONVERSION_PROB, then mother's genotype is CC
        motherHaplo.A = (1-FEMALE_CONVERSION_PROB)*0.5 + FEMALE_CONVERSION_PROB;
        motherHaplo.a = (1-FEMALE_CONVERSION_PROB)*0.5;
    else 
        % homing is not possible
        motherHaplo.A = motherHaplo.A + count(momGenotype,"A")/2;
        motherHaplo.a = motherHaplo.a + count(momGenotype,"a")/2;
    end    

    % fill in each column 
    % AA
    genoArray(i,1) = (fatherHaplo.A)*(motherHaplo.A);
    % Aa
    genoArray(i,2) = (fatherHaplo.A)*(motherHaplo.a) + (fatherHaplo.a)*(motherHaplo.A);
    % aa
    genoArray(i,3) = (fatherHaplo.a)*(motherHaplo.a);

end

% append a column for death probability
genoArray = [genoArray, 1-sum(genoArray,2)];
% correct numerical error
genoArray(abs(genoArray) < 10^(-4)) = 0; 


%% return
breedingTables = struct();
breedingTables.genoArray = genoArray;

end








