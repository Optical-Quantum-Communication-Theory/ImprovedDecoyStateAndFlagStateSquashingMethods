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

function channelModel = pmSixStateWCPChannel_constr_nodecoy(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["ed","pz","pd","eta","etad","mu1","active","fullstat"];
        
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
    
    decoys = [mu1];
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
    
    %save("BipartiteExps.mat","bipartiteExpectationsWCP","ObsFlag")
    %disp(ObsFlag)

    %Perform decoy analysis 
    %Define empty lists to store upper and lower bounds
    L1=size(bipartiteExpectationsWCP,1);
    L2=size(bipartiteExpectationsWCP,2);
    bipartiteExpectationsL=zeros(L1,L2);
    bipartiteExpectationsU=zeros(L1,L2);

    %Run decoy analysis    
    %Use simple LP for decoy
    [bipartiteExpectationsL,bipartiteExpectationsU] = decoyAnalysis(decoys,bipartiteExpectationsWCP);


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

%helper function that performs decoy state analysis
function [Y1L,Y1U] = decoyAnalysis(decoys,decoy_expecfull)
    %CVX solver
    

    m = size(decoy_expecfull,1);
    n = size(decoy_expecfull,2);

    Y1L = zeros(m,n);
    Y1U = zeros(m,n);

    n_photon=10;
    n_decoy=size(decoys,2);

    Poisson=@(mu,n) exp(-mu).*mu.^n./factorial(n);

    decoy_tolerance=1e-12;

    %Pick n=0,1,...-photon components
    M = {1,n_photon+1};
    for i = 1:n_photon+1
        Mi = zeros(m, m*(n_photon + 1));
        indr = [1:m]';
        indc = [i:n_photon+1:m*(n_photon + 1)]';
        indx = sub2ind(size(Mi),indr,indc);
        Mi(indx) = 1;
        M{i}=Mi;
    end
    
    %solve for upper bound
    for i=0:m-1
        for j=1:n
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            try
                cvx_begin quiet
                cvx_precision medium
                cvx_solver Mosek
                    variable Y(m*(n_photon+1),n)
                    maximize Objrow*Y*Objcolumn
                    subject to
                        Y >= 0;
                        Y <= 1;
                                                
                        Y0 = M{1}*Y;
                        
                        %e_0 =1/2
                        Y0(1,2) + Y0(2,1) == 1/2*(Y0(1,1)+Y0(1,2)+Y0(2,1)+Y0(2,2));
                        Y0(3,4) + Y0(4,3) == 1/2*(Y0(3,3)+Y0(3,4)+Y0(4,3)+Y0(4,4));
                        Y0(5,6) + Y0(6,5) == 1/2*(Y0(5,5)+Y0(5,6)+Y0(6,5)+Y0(6,6));

                        for k = 1:n_decoy
                            pmu = Poisson(decoys(k),0:1:n_photon);
                            Pmu = kron(eye(m),pmu);
                            Pmu*Y >= decoy_expecfull(:,:,k) - ones(m,n) * (1-sum(pmu)) - decoy_tolerance;
                            Pmu*Y <= decoy_expecfull(:,:,k) + decoy_tolerance;                            
                        end

                        for l_ph=1:n_photon+1
                            Yn = M{l_ph}*Y;                            
                            sum(sum(Yn)) == 4;                            
                            for l=1:m
                                sum(Yn(l,:)) == 1;
                            end
                        end
                cvx_end
                if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            catch 
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            Y1U(i+1,j) = Y(i*(n_photon+1)+2,j);            
        end
    end
        
    %solve for lower bound
    for i=0:m-1
        for j=1:n
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            try
                cvx_begin quiet
                cvx_precision medium
                cvx_solver Mosek
                    variable Y(m*(n_photon+1),n)
                    minimize Objrow*Y*Objcolumn
                    subject to
                        Y >= 0;
                        Y <= 1;

                        Y0 = M{1}*Y;
                        
                        %e_0 =1/2
                        Y0(1,2) + Y0(2,1) == 1/2*(Y0(1,1)+Y0(1,2)+Y0(2,1)+Y0(2,2));
                        Y0(3,4) + Y0(4,3) == 1/2*(Y0(3,3)+Y0(3,4)+Y0(4,3)+Y0(4,4));
                        Y0(5,6) + Y0(6,5) == 1/2*(Y0(5,5)+Y0(5,6)+Y0(6,5)+Y0(6,6));

                        for k = 1:n_decoy
                            pmu = Poisson(decoys(k),0:1:n_photon);
                            Pmu = kron(eye(m),pmu);
                            Pmu*Y >= decoy_expecfull(:,:,k) - ones(m,n) * (1-sum(pmu)) - decoy_tolerance;
                            Pmu*Y <= decoy_expecfull(:,:,k) + decoy_tolerance;                            
                        end

                        for l_ph=1:n_photon+1
                            Yn = M{l_ph}*Y;
                            sum(sum(Yn)) == 4;  
                            for l=1:m
                                sum(Yn(l,:)) == 1;
                            end                           
                        end
                cvx_end
                if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            catch 
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            Y1L(i+1,j) = Y(i*(n_photon+1)+2,j);
        end
    end
end

%helper function that simulates the channel for WCP sources
function expectations=coherentSourceChannel(active,mu,eta,etad,pd,ed,px)
    expectations=zeros(4,16);
    pz=1-px;
    t=eta*etad; %total transmittance
    Poisson=@(mu,n) exp(-mu)*mu^n/factorial(n);

    theta=asin(sqrt(ed));
    PL=sin(pi/4-theta)^2;
    PU=cos(pi/4-theta)^2;
    
    mapping_passive=[[pz*(1-ed),pz*ed,px*PU,px*PL];[pz*ed,pz*(1-ed),px*PL,px*PU];[pz*PL,pz*PU,px*(1-ed),px*ed];[pz*PU,pz*PL,px*ed,px*(1-ed)]];
    mapping_active=[[1-ed,ed,PU,PL];[ed,1-ed,PL,PU];[PL,PU,1-ed,ed];[PU,PL,ed,1-ed]];
    
    for input=1:4
        %iterating over each input state
        
        for output=1:16
            %iterating over each pattern
            a=index1to4(output-1); %detector event, 4 elements corresponding to [H,V,D,A]
            
            Ppattern=1;
            if(active==0)
                %passive basis choice
                for k=1:4
                    %iterating over each detector
                    Pclick=1-Poisson(mu*t*mapping_passive(input,k),0); 
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    if(a(k)==1)
                        Ppattern=Ppattern*Pclick;
                    elseif(a(k)==0)
                        Ppattern=Ppattern*(1-Pclick);
                    end
                end
            else
                %active basis choice
                PpatternZ=1;
                prob_activeZ=[1,1,0,0];
                for k=1:4
                    %iterating over each detector
                    Pclick=(1-Poisson(mu*t*mapping_active(input,k),0));
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    Pclick=prob_activeZ(k)*Pclick; %effect of basis choice (active)
                    if(a(k)==1)
                        PpatternZ=PpatternZ*Pclick;
                    elseif(a(k)==0)
                        PpatternZ=PpatternZ*(1-Pclick);
                    end
                end
                PpatternX=1;
                prob_activeX=[0,0,1,1];
                for k=1:4
                    %iterating over each detector
                    Pclick=(1-Poisson(mu*t*mapping_active(input,k),0));
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    Pclick=prob_activeX(k)*Pclick; %effect of basis choice (active)
                    if(a(k)==1)
                        PpatternX=PpatternX*Pclick;
                    elseif(a(k)==0)
                        PpatternX=PpatternX*(1-Pclick);
                    end
                end
                Ppattern=px*PpatternX+pz*PpatternZ;
            end
            expectations(input,output)=Ppattern;
        end
    end
end

%takes in an index of 1 and convert to array of 4
function a=index1to4(index)
    a(1) = floor(index/8);
    index = mod(index,8);
    a(2) = floor(index/4);
    index = mod(index,4);
    a(3) = floor(index/2);
    index = mod(index,2);
    a(4) = index;
end