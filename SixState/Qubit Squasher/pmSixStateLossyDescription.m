%% FUNCTION NAME: SixStateLossyDescription
% Channel description for the 6 protocol
%%

function protocolDescription = pmSixStateLossyDescription(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pzA","pxA", "pzB", "pxB", 'fullstat'];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this description file
    %can be automatically filled in by calling addObservables(x) or addObservables(x,'mask',maskValue)
    observables = {};
    obsMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description begin %%%%%%%%%%%%%%%%%%%%%%%%%
    %%fill in calculation process here

    dimA = 2;
    dimB = 3;
    
    numBases = 3; % x,y and z, sent by Alice
    
    %probability of Y-basis for Alice and Bob
    pyA = 1-pzA-pxA;
    pyB = 1-pzB-pxB;
    
    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, C the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
    krausOpZ = kron(kron( sqrt(pzA)*[[1,0];[0,1]], sqrt(pzB) * diag([1,1,0])), [1;0;0]); % for Z basis
    
    krausOpX = kron(kron( sqrt(pxA/2)*[[1,1];[1,-1]] , sqrt(pxB) * diag([1,1,0])), [0;1;0]); % for X basis
    
    krausOpY = kron(kron( sqrt(pyA/2)*[[1,1i];[1,-1i]], sqrt(pyB) * diag([1,1,0])), [0;0;1]); % for Y basis
    
    krausOp = {krausOpZ, krausOpX, krausOpY};
    

    %Check if Kraus operators sum to <= 1
    krausSum = 0;
    for i = 1:length(krausOp)
        krausSum = krausSum + krausOp{i}'*krausOp{i};
    end
    %disp(krausSum)

    % components for the pinching Z map
    % key projections (i.e. pinching channel)
    keyProj1 = kron(diag([1,0]), eye(dimB*numBases));
    keyProj2 = kron(diag([0,1]), eye(dimB*numBases));
    keyMap={keyProj1, keyProj2};
    

    %% Observables
    % Guarantee Alice's state is unchanged
    basisA = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basisA)
        addObservables(kron(basisA{iBasisElm}, eye(dimB)), 'mask', 0);
    end
    
    % normalization
    addObservables(eye(dimA*dimB), 'mask', 0);
    
    
    %Constructing the constraints in terms of bipartite POVMs
    %Relevant states
    h = [1;0];              v = [0;1];
    d = (h+v)/sqrt(2);      a = (h-v)/sqrt(2);
    r = (h+1i*v)/sqrt(2);   l = (h-1i*v)/sqrt(2);
    
    %Alice's POVMs
    basicAlicePOVMs = {pzA*[[1,0];[0,0]], pzA*[[0,0];[0,1]], pxA/2*[[1,1];[1,1]], pxA/2*[[1,-1];[-1,1]], pyA/2*[[1,1i];[-1i,1]], pyA/2*[[1,-1i];[1i,1]]};
    
    %Number of Alice's POVM elements
    dimPA = size(basicAlicePOVMs,2);
    
    %Check if Alice's POVMs sum to identity
    AlicePOVMSum = 0;
    
    for i = 1:length(basicAlicePOVMs)
        AlicePOVMSum = AlicePOVMSum + basicAlicePOVMs{i};
    end
    %disp(AlicePOVMSum)
    

    %Bob's POVMs
    basicBobPOVMs = {blkdiag(pzB*(h*(h')),0), blkdiag(pzB*(v*(v')),0), blkdiag(pxB*(d*(d')),0), blkdiag(pxB*(a*(a')),0), blkdiag(pyB*(r*(r')),0), blkdiag(pyB*(l*(l')),0), diag([0,0,1])};
    
    %Number of Bob's POVM elements
    dimPB = size(basicBobPOVMs,2);
    
    %Check if Alice's POVMs sum to identity
    BobPOVMSum = 0;
    
    for i = 1:length(basicBobPOVMs)
        BobPOVMSum = BobPOVMSum + basicBobPOVMs{i};
    end
    %disp(BobPOVMSum)
    
    %Full set of bipartite POVMS
    bipartitePOVMs = cell(dimPA*dimPB,1);
    for i = 1:dimPA
        for j = 1:dimPB
            bipartitePOVMs{dimPB*(i-1)+(j-1)+1} = kron(basicAlicePOVMs{i},basicBobPOVMs{j});
        end
    end
    
    
    if(fullstat==1)
        addObservables(bipartitePOVMs,'mask',1);
    else
        %QBER and Gain statistics
        %Add gain statistics
        disp("Add QBER and gain statistics");
    end
    
    %Check if observables are hermitian
    isherm = zeros(length(observables),1);
    for index=1:length(observables)
        isherm(index) = ishermitian(observables{index});
    end
    %disp(isherm);

    fprintf("number of observables: %d\n", length(observables));
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%

    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB];
end