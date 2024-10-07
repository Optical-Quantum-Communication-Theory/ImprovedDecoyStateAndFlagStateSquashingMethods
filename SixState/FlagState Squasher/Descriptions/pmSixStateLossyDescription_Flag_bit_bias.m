%% FUNCTION NAME: SixStateLossyDescription
% Channel description for the 6 protocol
%%

function protocolDescription = pmSixStateLossyDescription_Flag_bit_bias(names,p)

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
    
    dimA = 6;
    dimB = 14;
    
    numBases = 3; % x,y and z, sent by Alice
    
    %% Kraus operators for G and Z map
    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, C the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
    
    % Z basis
    % Sum of Bob's POVM elements in Z basis
    sumBobz = zeros(1,14);
    sumBobz(1:2) = sqrt(pzB);
    sumBobz(5:6) = 1;
    
    krausOpZ = kron(kron( [[1,0,0,0,0,0];[0,1,0,0,0,0]] , diag(sumBobz)), [1;0;0]); % for Z basis
    
    % X basis
    % Sum of Bob's POVM elements in X basis
    sumBobx = zeros(1,14);
    sumBobx(1:2) = sqrt(pxB);
    sumBobx(7:8) = 1;
    
    krausOpX = kron(kron( [[0,0,1,0,0,0];[0,0,0,1,0,0]] , diag(sumBobx)), [0;1;0]); % for X basis
    
    % Y basis
    % Sum of Bob's POVM elements in Y basis
    sumBoby = zeros(1,14);
    sumBoby(1:2) = sqrt(1-pxB-pzB);
    sumBoby(9:10) = 1;
    
    krausOpY = kron(kron( [[0,0,0,0,1,0];[0,0,0,0,0,1]] , diag(sumBoby)), [0;0;1]); % for Y basis
    
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
    basicAlicePOVMs = {diag([1,0,0,0,0,0]),diag([0,1,0,0,0,0]),diag([0,0,1,0,0,0]),diag([0,0,0,1,0,0]),diag([0,0,0,0,1,0]),diag([0,0,0,0,0,1])};
    
    %Bob's POVMs
    %POVMs in order: 
    %single clicks: H, V, +, -, R, L
    %double clicks: HV, +-, RL
    %cross clicks
    %diag([0,1,0,0,0,0,0,0,0,0,0,0]) = vac, flag vac ,flag H, flag V, flag +, flag -,
    %flag R, flag L, flag double-cLick HV, flag double-cLick +-, flag double-cLick RL, 
    %flag cross-cLick

    %Order of POVM elements H, V, +, -, R, L, DC-HV, DC-+-, DC-RL, CC, vac
    basicBobPOVMs = {blkdiag(pzB*(h*(h')), diag([0,0,1,0,0,0,0,0,0,0,0,0])), ...
        blkdiag(pzB*(v*(v')), diag([0,0,0,1,0,0,0,0,0,0,0,0])), ...
        blkdiag(pxB*(d*(d')), diag([0,0,0,0,1,0,0,0,0,0,0,0])), ...
        blkdiag(pxB*(a*(a')), diag([0,0,0,0,0,1,0,0,0,0,0,0])), ...
        blkdiag((1-pxB-pzB)*(r*(r')), diag([0,0,0,0,0,0,1,0,0,0,0,0])), ...
        blkdiag((1-pxB-pzB)*(l*(l')), diag([0,0,0,0,0,0,0,1,0,0,0,0])), ... 
        diag([0,0,0,0,0,0,0,0,0,0,1,0,0,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,0,1,0,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,0,0,1,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,0,0,0,1]), ...
        diag([0,0,1,1,0,0,0,0,0,0,0,0,0,0])};
    
    dimPB = size(basicBobPOVMs,2);
    
    BobPOVMSum = 0;
    
    for i = 1:length(basicBobPOVMs)
        BobPOVMSum = BobPOVMSum + basicBobPOVMs{i};
    end
    
    %disp(BobPOVMSum)
    
    %Full set of bipartite POVMS
    bipartitePOVMs = cell(dimA*dimPB,1);
    for i = 1:dimA
        for j = 1:dimPB
            bipartitePOVMs{dimPB*(i-1)+(j-1)+1} = kron(basicAlicePOVMs{i},basicBobPOVMs{j});
        end
    end

    disp(numel(bipartitePOVMs))
    
    if(fullstat==1)
        addObservables(bipartitePOVMs,'mask',1);
    else
        %QBER and Gain statistics
        %Add gain statistics
        disp("Add QBER and gain statistics");
    end

    %Flag-state operators
    % Bob's projector onto <=1-photon subspace
    Proj_sub = blkdiag(eye(3),zeros(11));
    
    %Full set of bipartite Flag-POVMs
    bipartiteFlagOps = cell(dimA,1);
    for i = 1:dimA
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