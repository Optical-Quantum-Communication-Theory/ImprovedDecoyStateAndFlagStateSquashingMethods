%helper function that performs decoy state analysis
function [Y1L,Y1U] = decoyAnalysis_vec(decoys,decoy_expecfull)
    m = size(decoy_expecfull,1);
    n = size(decoy_expecfull,2);

    Y1L = zeros(m,n);
    Y1U = zeros(m,n);

    n_photon=10;
    n_decoy=size(decoys,2);

    Poisson=@(mu,n) exp(-mu).*mu.^n./factorial(n);

    decoy_tolerance=1e-12;
    
    %solve for upper bound
    disp("Solving for upper bounds")
    for i=0:m-1
        for j=1:n
            fprintf("\n Upper bounds solving entry %d of %d ",n*i+(j-1)+1,m*n)
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            try
                cvx_begin quiet
                cvx_solver Mosek
                cvx_precision default
                    variable Y(m*(n_photon+1),n)
                    maximize Objrow*Y*Objcolumn
                    subject to
                        Y >= 0;
                        Y <= 1;
                        for k = 1:n_decoy
                            pmu = Poisson(decoys(k),0:1:n_photon);
                            Pmu = kron(eye(m),pmu);
                            Pmu*Y >= decoy_expecfull(:,:,k) - ones(m,n) * (1-sum(pmu)) - decoy_tolerance;
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
            Y1U(i+1,j) = Y(i*(n_photon+1)+2,j);
        end
    end
    
    %solve for lower bound
    disp("Solving for lower bounds")
    for i=0:m-1
        for j=1:n
            fprintf("\n Lower bounds solving entry %d of %d ",n*i+(j-1)+1,m*n)

            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            try
                cvx_begin quiet
                cvx_solver Mosek
                cvx_precision default
                    variable Y(m*(n_photon+1),n)
                    minimize Objrow*Y*Objcolumn
                    subject to
                        Y >= 0;
                        Y <= 1;
                        for k = 1:n_decoy
                            pmu = Poisson(decoys(k),0:1:n_photon);
                            Pmu = kron(eye(m),pmu);
                            Pmu*Y >= decoy_expecfull(:,:,k) - ones(m,n) * (1-sum(pmu)) - decoy_tolerance;
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
            Y1L(i+1,j) = Y(i*(n_photon+1)+2,j);
        end
    end
end