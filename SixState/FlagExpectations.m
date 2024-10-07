function [bipartiteExpectationsWCP, MFlag] = FlagExpectations(signalIntensity,eta,etad,pd,ed,pz,px,py)
    %Models a lossy channel with misalignment.
    %taken from the channel model. This gives Bob's epxectations for
    %H,V,+,-,R,L,Vac conditioned on each choice of polarization Alice could
    %have used from H,V,+,-,R,L. (multiply by Alice's probabilities on the left
    %as a diagonal matrix to remove conditoning).
    
    %Bob's probabilities for each basis
    bProb = [pz, px, py];
    
    %Raw expectations
    rawExpectations = coherentSourceChannel2(signalIntensity,eta,etad,pd,ed,bProb);
    
    %Mapping to simplify raw expectations
    mapping = zeros(64,11);
    
    %detector order: | H | V | + | - | R | L | DC-HV | DC-+- | DC-RL | CC | vac |
    
    %the vast majority of click patterns are cross clicks which we map to the
    %cross click column
    mapping(:,10) = 1;
    
    %vacuum to vacuum
    mapping = quickMap(mapping,[1,1,1,1,1,1],[0,0,0,0,0,0,0,0,0,0,1]);
    
    %single clicks. They map back to themselves
    mapping = quickMap(mapping,[2,1,1,1,1,1],[1,0,0,0,0,0,0,0,0,0,0]); %H
    mapping = quickMap(mapping,[1,2,1,1,1,1],[0,1,0,0,0,0,0,0,0,0,0]); %V
    mapping = quickMap(mapping,[1,1,2,1,1,1],[0,0,1,0,0,0,0,0,0,0,0]); %+
    mapping = quickMap(mapping,[1,1,1,2,1,1],[0,0,0,1,0,0,0,0,0,0,0]); %-
    mapping = quickMap(mapping,[1,1,1,1,2,1],[0,0,0,0,1,0,0,0,0,0,0]); %R
    mapping = quickMap(mapping,[1,1,1,1,1,2],[0,0,0,0,0,1,0,0,0,0,0]); %L
    
    %double clicks (in the same basis), mapped to their on column
    mapping = quickMap(mapping,[2,2,1,1,1,1],[0,0,0,0,0,0,1,0,0,0,0]); %Z (HV)
    mapping = quickMap(mapping,[1,1,2,2,1,1],[0,0,0,0,0,0,0,1,0,0,0]); %X (+-)
    mapping = quickMap(mapping,[1,1,1,1,2,2],[0,0,0,0,0,0,0,0,1,0,0]); %Y (RL)

    % Bipartite expectations
    bipartiteExpectationsWCP = rawExpectations*mapping;

    %Obeservation used for flag
    MFlag = 1 - sum(bipartiteExpectationsWCP(:,1:6),2) - bipartiteExpectationsWCP(:,11);
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
    
            %construct the equivalent coherent state Alice sent
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

function coherentVec = coherentSourceVec(basis,key,intensity)
    %Prepare a coherent state polarization in the given basis and signal state
    coherentVec = sqrt(intensity)*pauliBasis(basis,false)*zket(2,key);
end

function vec = coherentRotation(vec, theta, axisZXY)
    vec = rotMatrix(theta,axisZXY)*vec;
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

function mapping = quickMap(mapping,pattern,remaping)
    mapping(sub2indPlus([2,2,2,2,2,2],pattern),:) = remaping;
end
