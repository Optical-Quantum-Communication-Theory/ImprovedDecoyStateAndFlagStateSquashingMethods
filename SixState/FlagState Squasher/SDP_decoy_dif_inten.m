%helper function that performs decoy state analysis
function [Y1L,Y1U] = SDP_decoy_dif_inten(decoys,pz,px,py,decoy_expecfull,ObsFlag)  
    
    %Dimension of Bob
    dimB_SDP = 13;
    
    %Dimension of Alice
    dimA_SDP = 3;
    
    %cut-off for decoy
    n_photon=10;

    %number of intensities used
    n_decoy=size(decoys,2);
    n_signal_int = size(decoys{1},2);

    %m = #number of rows of observations
    %n = #number of colums of observations
    m = size(decoy_expecfull,1);
    n = size(decoy_expecfull,2);
    
    Y1L = zeros(m,n);
    Y1U = zeros(m,n);
    
    %Possonian distribution
    Poisson=@(mu,n) exp(-mu).*mu.^n./factorial(n);
    
    %Tolerances for decoy and Choi matrices
    decoy_tolerance = 1e-10;
    choi_tolerance = 1e-8;
    
    %Signal states
    h = [1;0];              v = [0;1];
    d = (h+v)/sqrt(2);      a = (h-v)/sqrt(2);
    r = (h+1i*v)/sqrt(2);   l = (h-1i*v)/sqrt(2);
    
    %Signal density matrices
    rhoH = blkdiag(h*h',[0]);
    rhoV = blkdiag(v*v',[0]);
    rhoP = blkdiag(d*d',[0]);
    rhoM = blkdiag(a*a',[0]);
    rhoR = blkdiag(r*r',[0]);
    rhoL = blkdiag(l*l',[0]);
    rhovac = diag([0,0,1]);
    
    %combine in a list
    rhox = {rhoH,rhoV,rhoP,rhoM,rhoR,rhoL};
    
    %Bob's POVM elements. This should not be here for future
    %implementations
    BobMeas = {blkdiag(pz*(h*(h')), diag([0,1,0,0,0,0,0,0,0,0,0])), ...
        blkdiag(pz*(v*(v')), diag([0,0,1,0,0,0,0,0,0,0,0])), ...
        blkdiag(px*(d*(d')), diag([0,0,0,1,0,0,0,0,0,0,0])), ...
        blkdiag(px*(a*(a')), diag([0,0,0,0,1,0,0,0,0,0,0])), ...
        blkdiag(py*(r*(r')), diag([0,0,0,0,0,1,0,0,0,0,0])), ...
        blkdiag(py*(l*(l')), diag([0,0,0,0,0,0,1,0,0,0,0])), ... 
        diag([0,0,0,0,0,0,0,0,0,1,0,0,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,1,0,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,0,1,0]), ...
        diag([0,0,0,0,0,0,0,0,0,0,0,0,1]), ...
        diag([0,0,1,0,0,0,0,0,0,0,0,0,0])};
    
    %Projection onto n<=N_B subspace
    Proj_sub = blkdiag(eye(3),zeros(10));
            
    %Matrices to pick n=0,1,...-photon components
    M = {1,n_photon+1};
    for i = 1:n_photon+1
        Mi = zeros(m, m*(n_photon + 1));
        indr = [1:m]';
        indc = [i:n_photon+1:m*(n_photon + 1)]';
        indx = sub2ind(size(Mi),indr,indc);
        Mi(indx) = 1;
        M{i}=Mi;
    end
    
     %solve for upper bound 1-photon
    disp("Upper bound single photons")
    
    for i=0:m-1
        for j=1:n
            fprintf("\n Upper bounds 1-photon solving SDP %d of %d ",n*i+(j-1)+1,m*n)

            %Select the correct rows and columns from Y for the objective
            %function
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;

            %Set up SDP for upper bounds
            % Y is the matrix representing all yields
            cvx_begin sdp quiet
                cvx_solver Mosek
                cvx_precision medium
                variable Y(m*(n_photon+1),n)
                variable J(dimA_SDP*dimB_SDP,dimA_SDP*dimB_SDP) hermitian semidefinite
                minimize -Objrow*Y*Objcolumn
                subject to
                    % 0 <= Y <=1
                    vec(Y) >= vec(zeros(m*(n_photon+1),n));
                    vec(Y) <= vec(ones(m*(n_photon+1),n));
                    
                    %0- and 1-photon component treated seperately with add.
                    %constraints
                    Y0 = M{1}*Y;
                    Y1 = M{2}*Y;

                    %Partial trace of Choi matrix = id
                    norm(PartialTrace(J,[1],[dimB_SDP,dimA_SDP])-eye(dimA_SDP)) <= choi_tolerance;
                                            
                    %Usual decoy bounds rewritten as matrix times vector
                    for k = 1:n_decoy
                        %total probability for given intensity to have
                        %<= n_photon photons
                        %create empty vector which is filled for each
                        %intensity again
                        pmu_tot = zeros(m,n);

                        %itereate over each subintensity
                        for kint = 1: n_signal_int
                            %select subintensity
                            intensity = decoys{k};
                            %calculate prob. distribution
                            pmu = Poisson(intensity(kint),0:1:n_photon);
                            %create list with all subintensities
                            pmu_list{kint} = pmu;
                            pmu_tot(kint,:) = sum(pmu);
                        end
                        %create matrix with vector of pmu on diagonal
                        Pmu = blkdiag(pmu_list{:});

                        %Create constraints on yields
                        vec(Pmu*Y) >= vec(decoy_expecfull(:,:,k) - (1-pmu_tot) - decoy_tolerance);
                        vec(Pmu*Y) <= vec(decoy_expecfull(:,:,k) + decoy_tolerance);                            
                    end

                    %Additional constraints for 0- and 1-photon components
                    %in terms of Choi matrices
                    for indexrow = 1:m
                        for indexcolumn = 1:n
                            norm(Y1(indexrow,indexcolumn) - trace(J*kron(BobMeas{indexcolumn},transpose(rhox{indexrow})))) <= choi_tolerance;
                            norm(Y0(indexrow,indexcolumn) - trace(J*kron(BobMeas{indexcolumn},transpose(rhovac)))) <= choi_tolerance;
                        end
                    end

                    %Upper and lower bounds for n<=N_B subspace
                    for indexrow =1:m
                        intensity = decoys{1}(indexrow);
                        %1-photon
                        real(trace(J*kron(Proj_sub,transpose(rhox{indexrow})))) >= (1- 1/(Poisson(intensity,1))*ObsFlag(indexrow,1)/(1-sum([pz,px,py].^2))) - choi_tolerance;                       
                        real(trace(J*kron(Proj_sub,transpose(rhox{indexrow})))) <= 1 + choi_tolerance;
                        %0-photon
                        real(trace(J*kron(Proj_sub,transpose(rhovac)))) >= (1- 1/(Poisson(intensity,0))*ObsFlag(indexrow,1)/(1-sum([pz,px,py].^2))) - choi_tolerance;
                        real(trace(J*kron(Proj_sub,transpose(rhovac)))) <= 1 + choi_tolerance;
                    end
            cvx_end
            disp(cvx_status)

            %Store upper bound
            Y1U(i+1,j) = Y(i*(n_photon+1)+2,j);            
        end
    end
        
    %solve for lower bound
    for i=0:m-1
        for j=1:n
            fprintf("\n Lower bounds 1-photon solving SDP %d of %d ",n*i+(j-1)+1,m*n)

            %Select the correct rows and columns from Y for the objective
            %function
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            %Set up SDP for lower bounds
            % Y is the matrix representing all yields
            cvx_begin sdp quiet
                cvx_solver Mosek
                cvx_precision medium
                variable Y(m*(n_photon+1),n)
                variable J(dimA_SDP*dimB_SDP,dimA_SDP*dimB_SDP) hermitian semidefinite
                minimize Objrow*Y*Objcolumn
                subject to
                    % 0 <= Y <=1
                    vec(Y) >= vec(zeros(m*(n_photon+1),n));
                    vec(Y) <= vec(ones(m*(n_photon+1),n));
                    
                    %0- and 1-photon component treated seperately with add.
                    %constraints
                    Y0 = M{1}*Y;
                    Y1 = M{2}*Y;

                    %Partial trace of Choi matrix = id
                    norm(PartialTrace(J,[1],[dimB_SDP,dimA_SDP])-eye(dimA_SDP)) <= choi_tolerance;
                                            
                    %Usual decoy bounds rewritten as matrix times vector
                    for k = 1:n_decoy
                        %total probability for given intensity to have
                        %<= n_photon photons
                        %create empty vector which is filled for each
                        %intensity again
                        pmu_tot = zeros(m,n);

                        %itereate over each subintensity
                        for kint = 1: n_signal_int
                            %select subintensity
                            intensity = decoys{k};
                            %calculate prob. distribution
                            pmu = Poisson(intensity(kint),0:1:n_photon);
                            %create list with all subintensities
                            pmu_list{kint} = pmu;
                            pmu_tot(kint,:) = sum(pmu);
                        end
                        %create matrix with vector of pmu on diagonal
                        Pmu = blkdiag(pmu_list{:});

                        %Create constraints on yields
                        vec(Pmu*Y) >= vec(decoy_expecfull(:,:,k) - (1-pmu_tot) - decoy_tolerance);
                        vec(Pmu*Y) <= vec(decoy_expecfull(:,:,k) + decoy_tolerance);                            
                    end

                    %Additional constraints for 0- and 1-photon components
                    %in terms of Choi matrices
                    for indexrow = 1:m
                        for indexcolumn = 1:n
                            norm(Y1(indexrow,indexcolumn) - trace(J*kron(BobMeas{indexcolumn},transpose(rhox{indexrow})))) <= choi_tolerance;
                            norm(Y0(indexrow,indexcolumn) - trace(J*kron(BobMeas{indexcolumn},transpose(rhovac)))) <= choi_tolerance;
                        end
                    end

                    %Upper and lower bounds for n<=N_B subspace
                    for indexrow =1:m
                        intensity = decoys{1}(indexrow);
                        %1-photon
                        trace(J*kron(Proj_sub,transpose(rhox{indexrow}))) >= (1- 1/(Poisson(intensity,1))*ObsFlag(indexrow,1)/(1-sum([pz,px,py].^2))) - choi_tolerance;                       
                        trace(J*kron(Proj_sub,transpose(rhox{indexrow}))) <= 1 + choi_tolerance;
                        %0-photon
                        trace(J*kron(Proj_sub,transpose(rhovac))) >= (1- 1/(Poisson(intensity,0))*ObsFlag(indexrow,1)/(1-sum([pz,px,py].^2))) - choi_tolerance;
                        trace(J*kron(Proj_sub,transpose(rhovac))) <= 1 + choi_tolerance;
                    end
            cvx_end
            disp(cvx_status)

            %Store lower bound
            Y1L(i+1,j) = Y(i*(n_photon+1)+2,j);
        end
    end
end