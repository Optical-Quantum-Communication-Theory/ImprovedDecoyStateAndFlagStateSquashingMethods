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

function channelModel = pmSixStateWCPChannel_Flag_bit_bias(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["ed","pzA","pxA","pyA", "pzB", "pxB","pyB","pd","eta","etad","mu1","mu2","mu3","fullstat"];
    %varNames=["ed","pzA","pxA","pyA", "pzB", "pxB","pyB","pd","eta","etad","mu1","mu2","fullstat"];

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
    decoys = {mu1,mu2,mu3};
    %decoys = {mu1,mu2};
    
    %Define signal intesnity as first decoy intensity
    signal_int = decoys{1};
    
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

    %Probability of sending signal given 1-photon was sent, i.e. p(x|1)
    probList_single_photon = [pzA/2; pzA/2; pxA/2; pxA/2; pyA/2; pyA/2];

    %Contraints on Alice's sent states, i.e. rho_A constraints
    % here with the reduced system rho_A is just the identity

    % rho_A constraints
    rhoA = zeros(dimA);
    %partial trace over flying qubit system to obtain local rhoA
    for jRow = 1 : dimA
        for kColumn = 1 : dimA
            rhoA(jRow,kColumn) = sqrt(probList_single_photon(jRow) * probList_single_photon(kColumn)) * signalStates{kColumn}' * signalStates{jRow};
        end
    end

    %add Alice's constraints to expectatons
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(rhoA * basis{iBasisElm}),'mask',0);
    end

    % Normalization
    addExpectations(1,'mask',0);

    %simulating the channel (already squashed), including observations for
    %flag-state squasher (ObsFlag)

    n_signal_int = size(decoys{1},2);

    for kdecoy =1:length(decoys)
        decoyint = decoys{kdecoy};
        for kint = 1:n_signal_int
            [obs, flag_obs(:,kdecoy)] = FlagExpectations(decoyint(kint),eta,etad,pd,ed,pzB,pxB,pyB);
            bipartiteExpectationsWCP(kint,:,kdecoy)= obs(kint,:);
            ObsFlag(kint,kdecoy) = flag_obs(kint,kdecoy);
        end
    end                                                 
    
    %Perform decoy analysis 
    %Define empty lists to store upper and lower bounds
    L1=size(bipartiteExpectationsWCP,1);
    L2=size(bipartiteExpectationsWCP,2);
    bipartiteExpectationsL=zeros(L1,L2);
    bipartiteExpectationsU=zeros(L1,L2);

    %Run decoy analysis
    %Use simple decoy methods
    [bipartiteExpectationsL,bipartiteExpectationsU] = decoyAnalysis_dif_inten(decoys,bipartiteExpectationsWCP);

    %Use improved decoy methods
    % [bipartiteExpectationsL,bipartiteExpectationsU] = SDP_decoy_dif_inten(decoys,pzB,pxB,pyB,bipartiteExpectationsWCP,ObsFlag);

    %Calculate weight inside subspace for flag-state squasher (still
    %conditional on Alice's state)
    flagExpecU = sum(bipartiteExpectationsU(:,7:10),2);

    %normalize by signal state probabilities
    bipartiteExpectationsL = diag(probList_single_photon)*bipartiteExpectationsL;
    bipartiteExpectationsU = diag(probList_single_photon)*bipartiteExpectationsU;

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
        bipartiteFlagExpectationsL(i) = probList_single_photon(i) * (1- flagExpecU(i,1)/(1-sum([pzB,pxB,pyB].^2))); %change ^2 to ^(NB+1) for higher NB
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
    bipartiteFlagExpectationsU = probList_single_photon; 
    addExpectations(bipartiteFlagExpectationsU,'mask',2);
    
    %disp(expectations)
    
    %Probability of sending signal was sent, i.e. p(x)
    %Adjusted such that p(x|1)=1/2 => p_0 and p_1 change for each basis
    
    %Z basis
    pzA0 = Poisson(signal_int(2),1)/(Poisson(signal_int(1),1)+Poisson(signal_int(2),1));
    pzA1 = 1 -pzA0;
    
    %X basis
    pxA0 = Poisson(signal_int(4),1)/(Poisson(signal_int(3),1)+Poisson(signal_int(4),1));
    pxA1 = 1 -pxA0;

    %Y basis
    pyA0 = Poisson(signal_int(6),1)/(Poisson(signal_int(5),1)+Poisson(signal_int(6),1));
    pyA1 = 1 -pyA0;

    probList_signal = [pzA*pzA0; pzA*pzA1; pxA*pxA0; pxA*pxA1; pyA*pyA0; pyA*pyA1];


    %QBER and Gain statistics for error correction
    %Signal intensity is first intensity setting and used only for key generation
    
    signal_simulation = diag(probList_signal)*bipartiteExpectationsWCP(:,:,1);
    
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
    P1 = Poisson(signal_int,1)*probList_signal;
    P0 = Poisson(signal_int,0)*probList_signal;

    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx, errory, errorz];
    channelModel.pSift = [gainx, gainy, gainz];
    channelModel.pSignal = P1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prob = Poisson(mu,n)
    prob = exp(-mu).*mu.^n/factorial(n);
end