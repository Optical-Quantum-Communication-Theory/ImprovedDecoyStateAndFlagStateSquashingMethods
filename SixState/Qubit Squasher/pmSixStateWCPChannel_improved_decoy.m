% %% FUNCTION NAME: pmBB84WCPChannel
% Realistic channel model for prepare-and-measure BB84 using WCP source and decoy states. 
% The transmittance eta, misalignment ed, and dark count pd are considered.
% The expectations correspond to a squashing model with five POVM outcomes (including photon loss).
%
% The flags active and fullstat can switch between active/passive
% detection, and coarse-grain/fine-grain statistics.
%
% Decoy state analysis is performed inside this channel model function.
%
% Additional information including masks (denoting which statistics are
% bounded by decoy state analysis and therefore uncertain) 
% and pSignal (denoting single photon probability) are included.
%%

function channelModel = pmSixStateWCPChannel_improved_decoy(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["ed","pzA","pxA", "pzB", "pxB","pd","eta","etad","mu1","mu2","mu3","active","fullstat"];
    % varNames=["ed","pzA","pxA", "pzB", "pxB","pd","eta","etad","mu1","mu2","active","fullstat"];

    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this channel model file
    %can be automatically filled in by calling addExpectations(x) or addExpectations(x,'mask',maskValue)
    expectations = [];
    expMask = [];

    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    decoys = [mu1,mu2,mu3];
    % decoys = [mu1,mu2];
    nDecoys = length(decoys);

    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    
    dimPB = 7;

    pyA = 1 - pzA - pxA;
    pyB = 1 - pzB - pxB;

    %Signal states
    h = [1;0];              v = [0;1];
    d = (h+v)/sqrt(2);      a = (h-v)/sqrt(2);
    r = (h+1i*v)/sqrt(2);   l = (h-1i*v)/sqrt(2);

    signalStates = {h, v, d, a, r, l};
    probList = [pzA/2; pzA/2; pxA/2; pxA/2; pyA/2; pyA/2];

    %Contrainst on Alice's sent states, i.e. rho_A constraints
    rhoA = 1/2*eye(dimA);
    
    %add Alice's constraints to expectatons
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(rhoA * basis{iBasisElm}),'mask',0);
    end

    % Normalization
    addExpectations(1,'mask',0);

    %simulating the channel already squashed
    for i=1:length(decoys)
        bipartiteExpectationsWCP(:,:,i)=squashedExpectations(decoys(i),eta,etad,pd,ed,pzB,pxB,"coherent");
    end

    %Perform decoy analysis 
    decoy_expectations = squeeze(bipartiteExpectationsWCP);
    [bipartiteExpectationsL,bipartiteExpectationsU] = SDP_decoy_SixState(decoys,pxB,pyB,pzB,decoy_expectations);
       
    
    %normalize by signal state probabilities
    bipartiteExpectationsL = diag(probList)*bipartiteExpectationsL;
    bipartiteExpectationsU = diag(probList)*bipartiteExpectationsU;

    disp(bipartiteExpectationsL)
    disp(bipartiteExpectationsU)

    %Convert expectations after decoy analysis into vectors
    bipartiteExpectationsL_1D = zeros(6*7,1);
    bipartiteExpectationsU_1D = zeros(6*7,1);
    for i = 1:6
        for j = 1:7
            bipartiteExpectationsL_1D(7*(i-1)+(j-1)+1) = bipartiteExpectationsL(i,j);
            bipartiteExpectationsU_1D(7*(i-1)+(j-1)+1) = bipartiteExpectationsU(i,j);
        end
    end

    if(fullstat==1)
        addExpectations(bipartiteExpectationsL_1D,'mask',1);
        addExpectations(bipartiteExpectationsU_1D,'mask',2);
    else
        %QBER and Gain statistics
        %Add gain statistics
        disp("Add QBER and gain statistics");
    end

    
    %QBER and Gain statistics for error correction
    %Signal intensity is first intensity setting and used only for key generation 
    signal_simulation = diag(probList)*bipartiteExpectationsWCP(:,:,1);

    gainz=sum(signal_simulation.*[1,1,0,0,0,0,0; 1,1,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0],'all');
    errorz=sum(signal_simulation.*[0,1,0,0,0,0,0; 1,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0],'all')/gainz;
    gainx=sum(signal_simulation.*[0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,1,1,0,0,0; 0,0,1,1,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0],'all');
    errorx=sum(signal_simulation.*[0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,1,0,0,0; 0,0,1,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0],'all')/gainx;
    gainy=sum(signal_simulation.*[0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,1,1,0; 0,0,0,0,1,1,0],'all');
    errory=sum(signal_simulation.*[0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,0,0; 0,0,0,0,0,1,0; 0,0,0,0,1,0,0],'all')/gainy;


    %signal state proportion
    P1 = decoys(1)*exp(-decoys(1));


    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx, errory, errorz];
    channelModel.pSift = [gainx, gainy, gainz];
    channelModel.pSignal = P1;
end