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

function channelModel = pmBB84WCPChannel_constr_nodecoy(protocolDescription,names,p)

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
    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    nDecoys = length(decoys);
    dimPB = 5;
    px = 1 - pz;
    
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];
    signalStates = {[1;0], [0;1], ketPlus, ketMinus};
    probList = [pz/2; pz/2; (1-pz)/2; (1-pz)/2];
    
    % rho_A constraints
    rhoA = zeros(dimA);
    %partial trace over flying qubit system to obtain local rhoA
    for jRow = 1 : dimA
        for kColumn = 1 : dimA
            rhoA(jRow,kColumn) = sqrt(probList(jRow) * probList(kColumn)) * signalStates{kColumn}' * signalStates{jRow};
        end
    end
    expectations = [];
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(rhoA * basis{iBasisElm}),'mask',0);
    end
    
    % Normalization
    addExpectations(1,'mask',0);
    
    
    %simulating the channel
    for i=1:length(decoys)
        rawExpectations(:,:,i)=coherentSourceChannel(active,decoys(i),eta,etad,pd,ed,px);
    end
    
    %convert 16-D pattern to 5-D POVM
    %Patterns:
    %column: [H,V,+,-,None]
    %row: Alice POVM
    MappingH=zeros(16,1);
    MappingH(9)=1; %1000
    MappingH(13)=0.5; %1100
    MappingV=zeros(16,1);
    MappingV(5)=1; %0100
    MappingV(13)=0.5; %1100
    
    MappingD=zeros(16,1);
    MappingD(3)=1; %0010
    MappingD(4)=0.5; %0011
    
    MappingA=zeros(16,1);
    MappingA(2)=1; %0001
    MappingA(4)=0.5; %0011
    
    MappingN=zeros(16,1);
    MappingN(1)=1; %0000
    MappingN(6)=1; %0101
    MappingN(7)=1; %0110
    MappingN(8)=1; %0111
    MappingN(10)=1; %1001
    MappingN(11)=1; %1010
    MappingN(12)=1; %1011
    MappingN(14)=1; %1101
    MappingN(15)=1; %1110
    MappingN(16)=1; %1111
    
    Mapping=[MappingH,MappingV,MappingD,MappingA,MappingN];
    
    %first bin into squashed model, then perform decoy analysis
    for i=1:length(decoys)
        bipartiteExpectationsWCP(:,:,i)=rawExpectations(:,:,i)*Mapping;
    end

    % save("observations_single.mat","bipartiteExpectationsWCP")
    
    decoy_expectations=squeeze(bipartiteExpectationsWCP);
    [bipartiteExpectationsL,bipartiteExpectationsU]=decoyAnalysis(decoys,decoy_expectations);

    disp(bipartiteExpectationsL)
    disp(bipartiteExpectationsU)

    %normalize by signal state probabilities
    bipartiteExpectationsL = diag(probList)*bipartiteExpectationsL;
    bipartiteExpectationsU = diag(probList)*bipartiteExpectationsU;
    
    bipartiteExpectationsL_1D = zeros(4*5,1);
    bipartiteExpectationsU_1D = zeros(4*5,1);
    for i = 1:4
        for j = 1:5
            bipartiteExpectationsL_1D(5*(i-1)+(j-1)+1) = bipartiteExpectationsL(i,j);
            bipartiteExpectationsU_1D(5*(i-1)+(j-1)+1) = bipartiteExpectationsU(i,j);
        end
    end
    
    if(fullstat==1)
        addExpectations(bipartiteExpectationsL_1D,'mask',1);
        addExpectations(bipartiteExpectationsU_1D,'mask',2);
    else
        %QBER and Gain statistics
        select=@(x,y)dimPB*(x-1)+(y-1)+1;
        coarseGrainExpectationsL = [bipartiteExpectationsL_1D(select(1,2));bipartiteExpectationsL_1D(select(2,1));bipartiteExpectationsL_1D(select(1,1));bipartiteExpectationsL_1D(select(2,2));...
                bipartiteExpectationsL_1D(select(3,4));bipartiteExpectationsL_1D(select(4,3));bipartiteExpectationsL_1D(select(3,3));bipartiteExpectationsL_1D(select(4,4))];
        coarseGrainExpectationsU = [bipartiteExpectationsU_1D(select(1,2));bipartiteExpectationsU_1D(select(2,1));bipartiteExpectationsU_1D(select(1,1));bipartiteExpectationsU_1D(select(2,2));...
                bipartiteExpectationsU_1D(select(3,4));bipartiteExpectationsU_1D(select(4,3));bipartiteExpectationsU_1D(select(3,3));bipartiteExpectationsU_1D(select(4,4))];
        %normalization
        normL = 1-sum(coarseGrainExpectationsU);
        normU = 1-sum(coarseGrainExpectationsL);
        
        addExpectations(coarseGrainExpectationsL,'mask',1);
        addExpectations(normL,'mask',1);
        addExpectations(coarseGrainExpectationsU,'mask',2);
        addExpectations(normU,'mask',2);
    end
    
        
    %QBER and Gain statistics for error correction
    signal_simulation = diag(probList)*(rawExpectations(:,:,1)*Mapping);
    gainz=sum(signal_simulation.*[1,1,0,0,0;1,1,0,0,0;0,0,0,0,0;0,0,0,0,0],'all');
    errorz=sum(signal_simulation.*[0,1,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,0],'all')/gainz;
    gainx=sum(signal_simulation.*[0,0,0,0,0;0,0,0,0,0;0,0,1,1,0;0,0,1,1,0],'all');
    errorx=sum(signal_simulation.*[0,0,0,0,0;0,0,0,0,0;0,0,0,1,0;0,0,1,0,0],'all')/gainx;
    
    %signal state proportion
    P1 = decoys(1)*exp(-decoys(1));

    % save("signal_simul_single.mat","signal_simulation")
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx,errorz];
    channelModel.pSift = [gainx,gainz];
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
                            
                            Yn(1,2) + Yn(2,1) == 0;
                            Yn(3,4) + Yn(4,3) == 0;
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
                            
                            Yn(1,2) + Yn(2,1) == 0;
                            Yn(3,4) + Yn(4,3) == 0;
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