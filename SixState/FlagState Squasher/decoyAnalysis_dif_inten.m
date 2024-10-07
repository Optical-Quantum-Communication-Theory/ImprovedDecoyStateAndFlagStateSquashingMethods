%helper function that performs decoy state analysis
function [Y1L,Y1U] = decoyAnalysis_dif_inten(decoys,decoy_expecfull)   
    %number of intensities used
    n_decoy=size(decoys,2);
    n_signal_int = size(decoys{1},2);

    %m = #number of rows of observations
    %n = #number of colums of observations
    m = size(decoy_expecfull,1);
    n = size(decoy_expecfull,2);
    
    %Empty vectors for upper and lower bounds
    Y1L = zeros(m,n);
    Y1U = zeros(m,n);
    
    %cut-off for decoy
    n_photon=10;
    
    %Possonian distribution
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
    
     %solve for upper bound 1-photon
    disp("Upper bound single photons")
    
    for i=0:m-1
        for j=1:n
            fprintf("\n Upper bounds 1-photon solving LP %d of %d ",n*i+(j-1)+1,m*n)

            %Select the correct rows and columns from Y for the objective
            %function
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            try
                %Set up LP for upper bounds
                % Y is the matrix representing all yields
                cvx_begin quiet
                    %CVX solver
                    cvx_solver Mosek
                    cvx_precision default
                    variable Y(m*(n_photon+1),n)
                    maximize Objrow*Y*Objcolumn
                    subject to
                        % 0 <= Y <=1       
                        Y >= 0;
                        Y <= 1;
                        
                        %0-photon yields
                        Y0 = M{1}*Y;
                        
                        %0-photon error rate e_0=1/2
                        % Y0(1,2) + Y0(2,1) == 1/2*(Y0(1,1)+Y0(1,2)+Y0(2,1)+Y0(2,2));
                        % Y0(3,4) + Y0(4,3) == 1/2*(Y0(3,3)+Y0(3,4)+Y0(4,3)+Y0(4,4));
                        
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
                            Pmu*Y >= decoy_expecfull(:,:,k) - (1-pmu_tot) - decoy_tolerance;
                            Pmu*Y <= decoy_expecfull(:,:,k) + decoy_tolerance;                            
                        end
                cvx_end
                disp(cvx_status)
            if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            catch 
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end

            %Store upper bound
            Y1U(i+1,j) = Y(i*(n_photon+1)+2,j);            
        end
    end
        
    %solve for lower bound
    for i=0:m-1
        for j=1:n
            fprintf("\n Lower bounds 1-photon solving LP %d of %d ",n*i+(j-1)+1,m*n)

            %Select the correct rows and columns from Y for the objective
            %function
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            try
                %Set up LP for upper bounds
                % Y is the matrix representing all yields
                cvx_begin quiet
                    cvx_solver Mosek
                    cvx_precision default
                    variable Y(m*(n_photon+1),n)
                    minimize Objrow*Y*Objcolumn
                    subject to

                        % 0 <= Y <=1
                        Y >= 0;
                        Y <= 1;
                        
                        %0-photon yields
                        Y0 = M{1}*Y;

                        %0-photon error rate e_0=1/2
                        % Y0(1,2) + Y0(2,1) == 1/2*(Y0(1,1)+Y0(1,2)+Y0(2,1)+Y0(2,2));
                        % Y0(3,4) + Y0(4,3) == 1/2*(Y0(3,3)+Y0(3,4)+Y0(4,3)+Y0(4,4));

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
                            Pmu*Y >= decoy_expecfull(:,:,k) - (1-pmu_tot) - decoy_tolerance;
                            Pmu*Y <= decoy_expecfull(:,:,k) + decoy_tolerance;                            
                        end
                cvx_end
                disp(cvx_status)
            if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            catch 
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end

            %Store lower bound
            Y1L(i+1,j) = Y(i*(n_photon+1)+2,j);
        end
    end
end