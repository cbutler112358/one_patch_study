% Function algorithmically producing breeding tables for a two-locus
% engineered underdominance tethered homing drive.
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


function [breedingTables] = genBT_TLU_tethered_homing(MALE_CONVERSION_PROB, ...
    FEMALE_CONVERSION_PROB)

%% produce breeding table for male zygotes
tmpArray = readcell('TLU-TH_breeding_table.xlsx');

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
    fatherHaplo.C = 0;
    fatherHaplo.c = 0;
    
    motherHaplo = fatherHaplo;    
    
    dadGenotype = string(tmpArray(1+i,1));
    momGenotype = string(tmpArray(1+i,2));

    % extract locus genotypes for each parent
    dadAGeno = extractBetween(dadGenotype,1,2);
    dadBGeno = extractBetween(dadGenotype,3,4);
    dadCGeno = extractBetween(dadGenotype,5,6);
    momAGeno = extractBetween(momGenotype,1,2);
    momBGeno = extractBetween(momGenotype,3,4);
    momCGeno = extractBetween(momGenotype,5,6);  


    % progeny are produced, check for homing in either parent
    homingDadBool = (count(dadBGeno,"B") > 0) & (dadCGeno == "Cc");
    homingMomBool = (count(momBGeno,"B") > 0) & (momCGeno == "Cc");

    if (homingDadBool)
        % with MALE_CONVERSION_PROB, then father's genotype is CC
        fatherHaplo.C = (1-MALE_CONVERSION_PROB)*0.5 + MALE_CONVERSION_PROB;
        fatherHaplo.c = (1-MALE_CONVERSION_PROB)*0.5;
        fatherHaplo.A = fatherHaplo.A + count(dadAGeno,"A")/2;
        fatherHaplo.a = fatherHaplo.a + count(dadAGeno,"a")/2;
        fatherHaplo.B = fatherHaplo.B + count(dadBGeno,"B")/2;
        fatherHaplo.b = fatherHaplo.b + count(dadBGeno,"b")/2;        
    else 
        % homing is not possible
        fatherHaplo.A = fatherHaplo.A + count(dadAGeno,"A")/2;
        fatherHaplo.a = fatherHaplo.a + count(dadAGeno,"a")/2;
        fatherHaplo.B = fatherHaplo.B + count(dadBGeno,"B")/2;
        fatherHaplo.b = fatherHaplo.b + count(dadBGeno,"b")/2;
        fatherHaplo.C = fatherHaplo.C + count(dadCGeno,"C")/2;
        fatherHaplo.c = fatherHaplo.c + count(dadCGeno,"c")/2;
    end
    
    if (homingMomBool)
        % with FEMALE_CONVERSION_PROB, then mother's genotype is CC
        motherHaplo.C = (1-FEMALE_CONVERSION_PROB)*0.5 + FEMALE_CONVERSION_PROB;
        motherHaplo.c = (1-FEMALE_CONVERSION_PROB)*0.5;
        motherHaplo.A = motherHaplo.A + count(momAGeno,"A")/2;
        motherHaplo.a = motherHaplo.a + count(momAGeno,"a")/2;
        motherHaplo.B = motherHaplo.B + count(momBGeno,"B")/2;
        motherHaplo.b = motherHaplo.b + count(momBGeno,"b")/2;        
    else 
        % homing is not possible
        motherHaplo.A = motherHaplo.A + count(momAGeno,"A")/2;
        motherHaplo.a = motherHaplo.a + count(momAGeno,"a")/2;
        motherHaplo.B = motherHaplo.B + count(momBGeno,"B")/2;
        motherHaplo.b = motherHaplo.b + count(momBGeno,"b")/2;
        motherHaplo.C = motherHaplo.C + count(momCGeno,"C")/2;
        motherHaplo.c = motherHaplo.c + count(momCGeno,"c")/2;
    end    

    % fill in each column 
    % AABBCC
    genoArray(i,1) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.C);
    % AABBCc
    genoArray(i,2) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.C);
    % AABBcc
    genoArray(i,3) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.c);
    % AABbCC
    genoArray(i,4) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.C) + ...
        (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.C);
    % AABbCc
    genoArray(i,5) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.C) + ...
        (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.C);    
    % AABbcc
    genoArray(i,6) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.c);
    % AAbbCC
    genoArray(i,7) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.C);
    % AAbbCc
    genoArray(i,8) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.C);
    % AAbbcc
    genoArray(i,9) = (fatherHaplo.A)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.c);
    % AaBBCC
    genoArray(i,10) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.C);
    % AaBBCc
    genoArray(i,11) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.C);
    % AaBBcc
    genoArray(i,12) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.c);
    % AaBbCC
    genoArray(i,13) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.C) + ...
        (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.C);
    % AaBbCc
    genoArray(i,14) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.C) + ...
        (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.C);
    % AaBbcc
    genoArray(i,15) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.c);    
    % AabbCC
    genoArray(i,16) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.C);
    % AabbCc
    genoArray(i,17) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.C);    
    % Aabbcc
    genoArray(i,18) = (fatherHaplo.A)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.A)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.c);
    % aaBBCC
    genoArray(i,19) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.C);    
    % aaBBCc
    genoArray(i,20) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.C);
    % aaBBcc
    genoArray(i,21) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.c);
    % aaBbCC
    genoArray(i,22) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.C);
    % aaBbCc
    genoArray(i,23) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.C) + ...
        (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.C);
    % aaBbcc
    genoArray(i,24) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.B)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.B)*(fatherHaplo.c)*(motherHaplo.c);
    % aabbCC
    genoArray(i,25) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.C);
    % aabbCc
    genoArray(i,26) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.C)*(motherHaplo.c) + ...
        (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.C);
    % aabbcc
    genoArray(i,27) = (fatherHaplo.a)*(motherHaplo.a)*(fatherHaplo.b)*(motherHaplo.b)*(fatherHaplo.c)*(motherHaplo.c);

end

% append a column for death probability
genoArray = [genoArray, 1-sum(genoArray,2)];
% correct numerical error
genoArray(abs(genoArray) < 10^(-4)) = 0; 


%% return
breedingTables = struct();
breedingTables.genoArray = genoArray;

end








