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

function channelModel = pmSixStateWCPChannel_Flag_improved(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["ed","pzA","pxA","pyA", "pzB", "pxB","pyB","pd","eta","etad","mu1","mu2","mu3","fullstat"];
    %varNames=["ed","pzA","pxA","pyA", "pzB", "pxB","pyB","pd","eta","etad","mu1","mu2","fullstat"];
    %varNames=["ed","pzA","pxA","pyA", "pzB", "pxB","pyB","pd","eta","etad","mu1","fullstat"];

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
    
    %Decoy state intensities and number
    decoys = [mu1,mu2,mu3];
    %decoys = [mu1,mu2];
    %decoys = [mu1];

    nDecoys = length(decoys);

    %Dimensions of Alice and Bob
    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    
    %Number of Alice's and Bob's POVM elements
    dimPA = 6;    
    dimPB = 11;

    %Signal states
    h = [1;0];              v = [0;1];
    d = (h+v)/sqrt(2);      a = (h-v)/sqrt(2);
    r = (h+1i*v)/sqrt(2);   l = (h-1i*v)/sqrt(2);

    signalStates = {h, v, d, a, r, l};
    probList = [pzA/2; pzA/2; pxA/2; pxA/2; pyA/2; pyA/2];

    %Contraints on Alice's sent states, i.e. rho_A constraints
    % here with the reduced system rho_A is just the identity

    rhoA = 1/2*eye(dimA);
    
    %add Alice's constraints to expectatons
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(rhoA * basis{iBasisElm}),'mask',0);
    end

    % Normalization
    addExpectations(1,'mask',0);

    %simulating the channel (already squashed), including observations for
    %flag-state squasher (ObsFlag)

    for i=1:length(decoys)
        [bipartiteExpectationsWCP(:,:,i), ObsFlag(:,i)]= FlagExpectations(decoys(i),eta,etad,pd,ed,pzB,pxB,pyB);
    end                                                  

    %Perform decoy analysis 
    %Define empty lists to store upper and lower bounds
    L1=size(bipartiteExpectationsWCP,1);
    L2=size(bipartiteExpectationsWCP,2);
    bipartiteExpectationsL=zeros(L1,L2);
    bipartiteExpectationsU=zeros(L1,L2);

    %Run decoy analysis
    %Use improved decoy methods
    
    decoy_expectations = squeeze(bipartiteExpectationsWCP);
    [bipartiteExpectationsL,bipartiteExpectationsU]=SDP_decoy_SixState_Flag(decoys,pzB,pxB,pyB,decoy_expectations,ObsFlag);  


    %Calculate weight inside subspace for flag-state squasher (still
    %conditional on Alice's state)
    flagExpecU = sum(bipartiteExpectationsU(:,7:10),2);
    
    %normalize by signal state probabilities
    bipartiteExpectationsL = diag(probList)*bipartiteExpectationsL;
    bipartiteExpectationsU = diag(probList)*bipartiteExpectationsU;

    %Convert expectations after decoy analysis into vectors
    bipartiteExpectationsL_1D = zeros(dimPA*dimPB,1);
    bipartiteExpectationsU_1D = zeros(dimPA*dimPB,1);
    for i = 1:dimPA
        for j = 1:dimPB
            bipartiteExpectationsL_1D(dimPB*(i-1)+(j-1)+1) = bipartiteExpectationsL(i,j);
            bipartiteExpectationsU_1D(dimPB*(i-1)+(j-1)+1) = bipartiteExpectationsU(i,j);
        end
    end
    
    %Add lower bounds for decoy expectations first. This needs to be done
    %for the flag-state squasher
    if(fullstat==1)
        addExpectations(bipartiteExpectationsL_1D,'mask',1);
    else
        %QBER and Gain statistics
        %Add gain statistics
        disp("Add QBER and gain statistics");
    end

    %Lower bound for flag-state expectations
    bipartiteFlagExpectationsL = zeros(dimPA,1);

    for i = 1:dimPA
        bipartiteFlagExpectationsL(i) = probList(i) * (1- flagExpecU(i,1)/(1-sum([pzB,pxB,pyB].^2))); %change ^2 to ^(NB+1) for higher NB
    end
    addExpectations(bipartiteFlagExpectationsL,'mask',1);


    %Add upper bounds for decoy expectations. This needs to be done
    %for the flag-state squasher
    if(fullstat==1)
        addExpectations(bipartiteExpectationsU_1D,'mask',2);
    else
        %QBER and Gain statistics
        %Add gain statistics
        disp("Add QBER and gain statistics");
    end
    
    %Trivial upper bound <= p(x) for flag-state expectations
    bipartiteFlagExpectationsU = probList; %ones(6,1);%
    addExpectations(bipartiteFlagExpectationsU,'mask',2);
    
    %disp(expectations)
    
    %QBER and Gain statistics for error correction
    %Signal intensity is first intensity setting and used only for key generation 
    signal_simulation = diag(probList)*bipartiteExpectationsWCP(:,:,1);
    
    %Gain and error in Z basis
    %Matrix for gain in Z basis
    Matgainz = zeros(dimPA,dimPB);
    Matgainz(1:2,1:2) = 1;

    %Matrix for error in Z basis
    Materrorz = zeros(dimPA,dimPB);
    Materrorz(1,2) = 1;
    Materrorz(2,1) = 1;

    gainz=sum(signal_simulation.*Matgainz,'all');
    errorz=sum(signal_simulation.*Materrorz,'all')/gainz;

    %Gain and error in X basis
    %Matrix for gain in X basis
    Matgainx = zeros(dimPA,dimPB);
    Matgainx(3:4,3:4) = 1;

    %Matrix for error in X basis
    Materrorx = zeros(dimPA,dimPB);
    Materrorx(3,4) = 1;
    Materrorx(4,3) = 1;

    gainx=sum(signal_simulation.*Matgainx,'all');
    errorx=sum(signal_simulation.*Materrorx,'all')/gainx;

    %Gain and error in Y basis
    %Matrix for gain in Y basis
    Matgainy = zeros(dimPA,dimPB);
    Matgainy(5:6,5:6) = 1;

    %Matrix for error in Y basis
    Materrory = zeros(dimPA,dimPB);
    Materrory(5,6) = 1;
    Materrory(6,5) = 1;

    gainy=sum(signal_simulation.*Matgainy,'all');
    errory=sum(signal_simulation.*Materrory,'all')/gainy;


    %signal state proportion
    P1 = decoys(1)*exp(-decoys(1));


    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx, errory, errorz];
    channelModel.pSift = [gainx, gainy, gainz];
    channelModel.pSignal = P1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%