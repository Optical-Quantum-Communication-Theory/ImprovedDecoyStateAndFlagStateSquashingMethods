function [bipartiteExpectationsWCP] = squashedExpectations(signalIntensity,eta,etad,pd,ed,pz,px,sourceType)
%Models a lossy channel with misalignment.
%taken from the channel model. This gives Bob's epxectations for
%H,V,+,-,R,L,Vac conditioned on each choice of polarization Alice could
%have used from H,V,+,-,R,L. (multiply by Alice's probabilities on the left
%as a diagonal matrix to remove conditoning).

% pz = 1/3; %THIS CANNOT BE CHANGED. PROVEN ONLY WORKS FOR 1/3
% bProb = [pz,(1-pz)/2,(1-pz)/2]; % SAME GOES FOR THIS
py = 1- pz - px;
bProb = [pz, px, py];

%% simulating the channel
switch sourceType
    case "coherent"
        rawExpectations = coherentSourceChannel2(signalIntensity,eta,etad,pd,ed,bProb);
    case "single photon"
        rawExpectations = singlePhotonSourceChannel(eta,etad,pd,ed);
    case "vacuum"
        rawExpectations = vacuumSourceChannel(pd);
    case "multi photon"
        %phase randomized coherent source conditioned on the multi photon
        %section
        rawExpectations = multiPhotonSourceChannel(signalIntensity,eta,etad,pd,ed,bProb);
    otherwise
        error("Not one of the listed source types")
end


%% convert 64-D pattern to 7-D POVM
mapping = constructSquashingMap(pz);

%first bin into squashed model, then perform decoy analysis
bipartiteExpectationsWCP = rawExpectations*mapping;
end


function expectations = multiPhotonSourceChannel(signalIntensity,eta,etad,pd,ed,bProb)
%phase randomized coherent source conditioned on the multi photon section.
%Thankfully this can be constructed by using the expectations for coherent,
%single photon and vacuum.

%probability that the source sends single and vacuum states.
p0 = exp(-signalIntensity);
p1 = signalIntensity.*exp(-signalIntensity);

%I'm so happy that everything is linear and that coherent and phase
%randomized states have the same measured quantities (for Bob's detectors).
expectations = (coherentSourceChannel2(signalIntensity,eta,etad,pd,ed,bProb)...
    -p1*singlePhotonSourceChannel(eta,etad,pd,ed)...
    -p0*vacuumSourceChannel(pd))/(1-p0-p1);

end

function expectations = vacuumSourceChannel(pd)
%Statistics before squashing given that Alice sends vacuum

% No matter what polarization Alice picks, she's still sending the same
% vacuum state. A pure loss channel (as well as missalignment) has no
% effect.
bobDetectorsGivenAlice = zeros(6,7);
bobDetectorsGivenAlice(:,7) = 1;

%Now we see how dark counts effect this. We can use the same dark count map
%used for the single photon source
darkMap = constructDarkCountMap(pd);
expectations = bobDetectorsGivenAlice*darkMap;

end

function expectations = singlePhotonSourceChannel(eta,etad,pd,ed)
%Statistics before squashing given that Alice sends single photons in H,V,+,-,R,L

%total signal transmitted
transmit = eta*etad;

%missalignment angle and matrix
%  using theta = asin(sqrt(ed))
rotMat = rotMatrix(asin(sqrt(ed)),[0,0,1]);

bobDetectorsGivenAlice = zeros(6,7);

%start by creating the signal Alice sends (as a qubit)
for basisA = 1:3
    for keyA = 1:2
        indexA = sub2indPlus([2,3],[keyA,basisA]);

        rhoIn = pauliBasis(basisA,false)*zket(2,keyA);

        %apply the missalignment rotation
        rhoIn = rotMat*rhoIn;

        for basisB =1:3
            probDetClicks = pauliBasis(basisB,false)'*rhoIn;
            %absolute value squared, not lost, and sent to the right detector
            probDetClicks = transmit/3*(probDetClicks.*conj(probDetClicks)); 

            for keyB = 1:2
                indexB = sub2indPlus([2,3],[keyB,basisB]);
                bobDetectorsGivenAlice(indexA,indexB) = probDetClicks(keyB);
            end
        end
    end
end

%now the last column is the probability that no detector fires (before dark
%counts). Becuase we have a single photon, the only way this can occur is
%if the photon is lost in transmition.
bobDetectorsGivenAlice(:,7) = 1- transmit;


%Now we have to work out dark counts and bring this from a single photon
%responce to the full detector click patterns, then apply to Bob's single
%photon detection patterns.
darkMap = constructDarkCountMap(pd);
expectations = bobDetectorsGivenAlice*darkMap;
end

function darkMap = constructDarkCountMap(pd)

%for the single photon source

% single and vacuum detector patterns
singleVacPats = {...
    [1,0,0,0,0,0],... %H
    [0,1,0,0,0,0],... %V
    [0,0,1,0,0,0],... %+
    [0,0,0,1,0,0],... %-
    [0,0,0,0,1,0],... %R
    [0,0,0,0,0,1],... %L
    [0,0,0,0,0,0],... %Vac
    };

numDetectors = 6;
dimPatterns=[2,2,2,2,2,2];
numPatterns = prod(dimPatterns);

darkMap = zeros(7,numPatterns);

for indexClickPat = 1:numPatterns
    clickPattern = ind2subPlus(dimPatterns,indexClickPat)-1;

    %Now we will compare the patterns to the ones Bob has
    for indexSingle = 1:numel(singleVacPats)
        darkCountPattern = clickPattern - singleVacPats{indexSingle};


        %make sure that the pattern can be achieved without forgeting
        %clicks.
        if any(darkCountPattern < 0,"all")
            transProb = 0;
        else
            %the mapping is physical and we have to determine how many
            %darcounts occure
            numDarkCounts = sum(darkCountPattern);
            numDets = sum(singleVacPats{indexSingle});

            transProb = pd^numDarkCounts*(1-pd)^(numDetectors-numDets-numDarkCounts);
        end

        %write the probability into the map
        darkMap(indexSingle,indexClickPat) = transProb;
    end
end
end

function expectations = coherentSourceChannel2(mu,eta,etad,pd,ed,bProb)
%construct the table  of the probability of any one of Bob's detectors
%clicking for each phase randomized state Alice sends to Bob.
bobClicksGivenAlice = constructBobClicksGivenAlice(mu,eta,etad,pd,ed,bProb);

%take that table and construct the probability of each click pattern for
%all the detectors together.
expectations = constructClickPatternProbs(bobClicksGivenAlice);
end


function bobClicksGivenAlice = constructBobClicksGivenAlice(mu,eta,etad,pd,ed,bProb)
%total transmittance of the channel
transmit = eta*etad;
intensity = mu*transmit;

%missalignment as an angle
theta = asin(sqrt(ed));

% Now we construct the table for the probability of Bob's detector for
% basisB and detector keyB, clicks for an input in basisA and key keyA,
% for the the transformations applied to the state Alice sent.
bobClicksGivenAlice = zeros(6,6);

for basisA =1:3
    for keyA = 1:2

        %construct the equivolent coherent state Alice sent
        %(This already has the loss built into the intensity)
        cohVec = coherentSourceVec(basisA,keyA,intensity);

        %misalignment rotation
        cohVec = coherentRotation(cohVec,theta,[0,0,1]);


        %convert to a linear index
        indexA = sub2indPlus([2,3],[keyA,basisA]); %This order groups the same basis together and alternates key bits
        for basisB = 1:3

            %rotate From basisB to Z
            cohVecTemp = pauliBasis(basisB,false)'*cohVec;

            %Now we account for how much intensity is sent to this detector
            %basis.
            cohVecTemp = sqrt(bProb(basisB))*cohVecTemp;


            for keyB = 1:2
                %convert to a linear index
                indexB = sub2indPlus([2,3],[keyB,basisB]);

                %probability the detector doesn't click
                probNoClick = exp(-cohVecTemp(keyB)*cohVecTemp(keyB)');

                % Now we determine the probability of clicking (and by
                % extension not clicking) while also concidering dark
                % counts.
                probClick = 1-probNoClick*(1-pd);

                bobClicksGivenAlice(indexA,indexB) = probClick;
            end
        end
    end
end
end


function expectations = constructClickPatternProbs(bobClicksGivenAlice)
%now We'll set up each click pattern for Bob's detectors given Alice's
%input

%set up the array that will hold all the input and output patterns
expectations = zeros(2*3,2^(2*3)); %6 by 64

dimPatterns=[2,2,2,2,2,2];
numPatterns = prod(dimPatterns);

%toggle between clickProb and 1-clickProb.
clickProbSwitch =@(click,clickProb) (click==0).*(1-clickProb) + (click~=0).*clickProb;

for basisA = 1:3
    for keyA = 1:2
        %convert to a linear index
        indexA = sub2indPlus([2,3],[keyA,basisA]);

        for patternIndex = 1:numPatterns
            %convert the pattern to a linear index
            patternVec = ind2subPlus(dimPatterns,patternIndex)-1; %-1 to make it a binary vector.
            expectations(indexA,patternIndex) = prod(clickProbSwitch(patternVec,bobClicksGivenAlice(indexA,:)));
        end
    end

end

end


function coherentVec = coherentSourceVec(basis,key,intensity)
%Prepare a coherent state polarization in the given basis and signal state
coherentVec = sqrt(intensity)*pauliBasis(basis,false)*zket(2,key);
end

function rotMat = rotMatrix(theta,axisZXY)
%used some of Scott's code as a base for this
%normalize the axis
axisZXY = axisZXY/norm(axisZXY);

PauliZ = [1,0;0,-1];
PauliX = [0,1;1,0];
PauliY = [0,-1i;1i,0];

rotMat = cos(theta/2)*eye(2)-1i*sin(theta/2)...
    *(axisZXY(1)*PauliZ+axisZXY(2)*PauliX+axisZXY(3)*PauliY);
end

function vec = coherentRotation(vec, theta, axisZXY)
vec = rotMatrix(theta,axisZXY)*vec;
end


function mapping = constructSquashingMap(pz)

%From Lars' writeup on the 6 state biased squasher
%chance to discard cross clicks
pas3 = (1-pz).^3/4 +pz.^3;
pDiscard = 1-pas3/(2*(1-pas3));

mapping = zeros(64,7);

%detector order: H V + - R L Discard

%the vast majority of click patterns are cross clicks which we map with
%equal probability to any of the outcomes. The easist way to then set up
%the entire mapping is to start by setting everything as a cross click,
%then overwriting anthying that needs to be changed.

mapping(:) = (1-pDiscard)/6;
mapping(:,7) = pDiscard; %many of the cross clicks are discarded

%Now we overwrite any entre that isn't a cross click.
%Remeber, matlab indexs starting at 1 :(
    function mapping = quickMap (mapping,pattern,remaping)
        mapping(sub2indPlus([2,2,2,2,2,2],pattern),:) = remaping;
    end
%vacuum to discard
mapping = quickMap(mapping,[1,1,1,1,1,1],[0,0,0,0,0,0,1]);

%single clicks. They map back to themselves
mapping = quickMap(mapping,[2,1,1,1,1,1],[1,0,0,0,0,0,0]); %H
mapping = quickMap(mapping,[1,2,1,1,1,1],[0,1,0,0,0,0,0]); %V
mapping = quickMap(mapping,[1,1,2,1,1,1],[0,0,1,0,0,0,0]); %+
mapping = quickMap(mapping,[1,1,1,2,1,1],[0,0,0,1,0,0,0]); %-
mapping = quickMap(mapping,[1,1,1,1,2,1],[0,0,0,0,1,0,0]); %R
mapping = quickMap(mapping,[1,1,1,1,1,2],[0,0,0,0,0,1,0]); %L

%double clicks (in the same basis). Equally distributed back to their basis
%outcomes
mapping = quickMap(mapping,[2,2,1,1,1,1],[1/2,1/2,0,0,0,0,0]); %Z (HV)
mapping = quickMap(mapping,[1,1,2,2,1,1],[0,0,1/2,1/2,0,0,0]); %X (+-)
mapping = quickMap(mapping,[1,1,1,1,2,2],[0,0,0,0,1/2,1/2,0]); %Y (RL)
end