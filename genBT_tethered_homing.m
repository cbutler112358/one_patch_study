% Function algorithmically producing breeding tables for a homing drive
% tethered to a biallelic locus. 
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


function [breedingTables] = genBT_tethered_homing(MALE_CONVERSION_PROB, ...
    FEMALE_CONVERSION_PROB)

%% produce breeding table for male zygotes
tmpArray = readcell('TH_breeding_table.xlsx');

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
    fatherHaplo.B = 0;
    fatherHaplo.b = 0;
    
    motherHaplo = fatherHaplo;    
    
    dadGenotype = string(tmpArray(1+i,1));
    momGenotype = string(tmpArray(1+i,2));

    % extract locus genotypes for each parent
    dadAGeno = extractBetween(dadGenotype,1,2);
    dadBGeno = extractBetween(dadGenotype,3,4);
    momAGeno = extractBetween(momGenotype,1,2);
    momBGeno = extractBetween(momGenotype,3,4);


    % progeny are produced, check for homing in either parent
    homingDadBool = (count(dadAGeno,"A") > 0) & (dadBGeno == "Bb");
    homingMomBool = (count(momAGeno,"A") > 0) & (momBGeno == "Bb");

    if (homingDadBool)
        % with MALE_CONVERSION_PROB, then father's genotype is CC
        fatherHaplo.B = (1-MALE_CONVERSION_PROB)*0.5 + MALE_CONVERSION_PROB;
        fatherHaplo.b = (1-MALE_CONVERSION_PROB)*0.5;
        fatherHaplo.A = fatherHaplo.A + count(dadAGeno,"A")/2;
        fatherHaplo.a = fatherHaplo.a + count(dadAGeno,"a")/2;
    else 
        % homing is not possible
        fatherHaplo.A = fatherHaplo.A + count(dadAGeno,"A")/2;
        fatherHaplo.a = fatherHaplo.a + count(dadAGeno,"a")/2;
        fatherHaplo.B = fatherHaplo.B + count(dadBGeno,"B")/2;
        fatherHaplo.b = fatherHaplo.b + count(dadBGeno,"b")/2;
    end
    
    if (homingMomBool)
        % with FEMALE_CONVERSION_PROB, then mother's genotype is CC
        motherHaplo.B = (1-FEMALE_CONVERSION_PROB)*0.5 + FEMALE_CONVERSION_PROB;
        motherHaplo.b = (1-FEMALE_CONVERSION_PROB)*0.5;
        motherHaplo.A = motherHaplo.A + count(momAGeno,"A")/2;
        motherHaplo.a = motherHaplo.a + count(momAGeno,"a")/2;
    else 
        % homing is not possible
        motherHaplo.A = motherHaplo.A + count(momAGeno,"A")/2;
        motherHaplo.a = motherHaplo.a + count(momAGeno,"a")/2;
        motherHaplo.B = motherHaplo.B + count(momBGeno,"B")/2;
        motherHaplo.b = motherHaplo.b + count(momBGeno,"b")/2;
    end    

    % fill in each column 
    % AABB
    genoArray(i,1) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B);
    % AABb
    genoArray(i,2) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b) + ...
        (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B);
    % AAbb
    genoArray(i,3) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b);
    % AaBB
    genoArray(i,4) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B);    
    % AaBb
    genoArray(i,5) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b) + ...
        (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B);
    % Aabb
    genoArray(i,6) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b);
    % aaBB
    genoArray(i,7) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B);    
    % aaBb
    genoArray(i,8) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b) + ...
        (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B);
    % aabb
    genoArray(i,9) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b);
    

end

% append a column for death probability
genoArray = [genoArray, 1-sum(genoArray,2)];
% correct numerical error
genoArray(abs(genoArray) < 10^(-4)) = 0; 


%% return
breedingTables = struct();
breedingTables.genoArray = genoArray;

end








