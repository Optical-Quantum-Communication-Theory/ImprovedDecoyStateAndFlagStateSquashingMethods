%% FUNCTION NAME: SixStateLossyDescription
% Channel description for the 6 protocol
%%

function protocolDescription = pmSixStateLossyDescription_Flag_Reduced(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pzA","pxA","pyA", "pzB", "pxB","pyB", 'fullstat'];
    
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
    dimB = 13;
    
    numBases = 3; % x,y and z, sent by Alice
    
    %% Kraus operators for G and Z map
    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, C the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
    
    % Z basis
    % Sum of Bob's POVM elements in Z basis
    sumBobz = zeros(1,dimB);
    sumBobz(1:2) = sqrt(pzB);
    sumBobz(4:5) = 1;

    krausOpZ = kron(kron( sqrt(pzA)*(kron([1;0],[[1,0]; [0,0]]) + kron([0;1],[[0,0]; [0,1]])), diag(sumBobz)), [1;0;0]);

    % X basis
    % Sum of Bob's POVM elements in X basis
    sumBobx = zeros(1,dimB);
    sumBobx(1:2) = sqrt(pxB);
    sumBobx(6:7) = 1;
    
    krausOpX = kron(kron( sqrt(pxA)/2*(kron([1;0],[[1,1]; [1,1]]) + kron([0;1],[[1,-1]; [-1,1]])) , diag(sumBobx)), [0;1;0]);

    % Y basis
    % Sum of Bob's POVM elements in Y basis
    sumBoby = zeros(1,dimB);
    sumBoby(1:2) = sqrt(pyB);
    sumBoby(8:9) = 1;
    
    krausOpY = kron(kron( sqrt(pyA)/2*(kron([1;0],[[1,1i]; [-1i,1]]) + kron([0;1],[[1,-1i]; [1i,1]]) ) , diag(sumBoby)), [0;0;1]);
    
    %Full Kraus ops
    krausOp = {krausOpZ, krausOpX, krausOpY};
    
    %Check if Kraus operators sum to <= 1
    krausSum = 0;
    for i = 1:length(krausOp)
        krausSum = krausSum + krausOp{i}'*krausOp{i};
    end
    
    %disp(krausSum)

    keyProj1 = kron(diag([1,0]), eye(dimA*dimB*numBases));
    keyProj2 = kron(diag([0,1]), eye(dimA*dimB*numBases));

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
    %POVMs in order: 
    %single clicks: H, V, +, -, R, L
    %double clicks: HV, +-, RL
    %cross clicks
    %diag([0,1,0,0,0,0,0,0,0,0,0]) = vac, flag H, flag V, flag +, flag -,
    %flag R, flag L, flag double-cLick HV, flag double-cLick +-, flag double-cLick RL, 
    %flag cross-cLick

    %Order of POVM elements H, V, +, -, R, L, DC-HV, DC-+-, DC-RL, CC, vac
    basicBobPOVMs = {blkdiag(pzB*(h*(h')), diag([0,1,0,0,0,0,0,0,0,0,0])), ...
        blkdiag(pzB*(v*(v')), diag([0,0,1,0,0,0,0,0,0,0,0])), ...
        blkdiag(pxB*(d*(d')), diag([0,0,0,1,0,0,0,0,0,0,0])), ...
        blkdiag(pxB*(a*(a')), diag([0,0,0,0,1,0,0,0,0,0,0])), ...
        blkdiag(pyB*(r*(r')), diag([0,0,0,0,0,1,0,0,0,0,0])), ...
        blkdiag(pyB*(l*(l')), diag([0,0,0,0,0,0,1,0,0,0,0])), ... 
        diag([0,0,0,0,0,0,0,0,0,1,0,0,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,1,0,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,0,1,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,0,0,1]), ...
        diag([0,0,1,0,0,0,0,0,0,0,0,0,0])};
    
    %Number of Bob's POVM elements
    dimPB = size(basicBobPOVMs,2);
    
    %Check if Bob's POVMs sum to identity
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

    %Flag-state operators
    % Bob's projector onto <=1-photon subspace
    Proj_sub = blkdiag(eye(3),zeros(10));
    
    %Full set of bipartite Flag-POVMs
    bipartiteFlagOps = cell(dimPA,1);
    for i = 1:dimPA
        bipartiteFlagOps{i} = kron(basicAlicePOVMs{i},Proj_sub);
    end

    addObservables(bipartiteFlagOps,'mask',1);
    
    fprintf("number of observables: %d\n", length(observables));
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%

    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB];
end