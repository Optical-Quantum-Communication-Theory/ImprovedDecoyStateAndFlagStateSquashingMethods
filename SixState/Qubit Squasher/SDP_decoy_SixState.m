%helper function that performs decoy state analysis
function [Y1L,Y1U] = SDP_decoy_SixState(decoys,pxB,pyB,pzB,decoy_expecfull)
       
    %cut-off for decoy
    n_photon=10;

    %number of intensities used
    n_decoy=size(decoys,2);

    %m = #number of rows of observations
    %n = #number of colums of observations
    m = size(decoy_expecfull,1);
    n = size(decoy_expecfull,2);
    
    Y1L = zeros(m,n);
    Y1U = zeros(m,n);
    
    %Possonian distribution
    Poisson=@(mu,n) exp(-mu).*mu.^n./factorial(n);
    
    %Tolerances for decoy and Choi matrices
    decoy_tolerance = 1e-12;
    choi_tolerance = 1e-12;
    
    %Signal states
    h = [1;0];              v = [0;1];
    d = (h+v)/sqrt(2);      a = (h-v)/sqrt(2);
    r = (h+1i*v)/sqrt(2);   l = (h-1i*v)/sqrt(2);
    
    %Signal density matrices
    rhoH = h*h';
    rhoV = v*v';
    rhoP = d*d';
    rhoM = a*a';
    rhoR = r*r';
    rhoL = l*l';
    
    %combine in a cell array
    rhox = {rhoH,rhoV,rhoP,rhoM,rhoR,rhoL};
    
    %vacuum state
    vac = [0;0;1];
    
    %Bob's POVM elements. This should not be here for future
    %implementations
    FH = pzB*blkdiag(rhoH,[0]);
    FV = pzB*blkdiag(rhoV,[0]);
    FP = pxB*blkdiag(rhoP,[0]);
    FM = pxB*blkdiag(rhoM,[0]);
    FR = pyB*blkdiag(rhoR,[0]);
    FL = pyB*blkdiag(rhoL,[0]);
    Fvac = vac*vac';
    
    %combine in a list
    Fy = {FH,FV,FP,FM,FR,FL,Fvac};
        
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

    %solve for upper bound
    disp("Solving for upper bounds")
    for i=0:m-1
        for j=1:n
            fprintf("\n Upper bounds solving SDP %d of %d ",n*i+(j-1)+1,m*n)
            
            %Select the correct rows and columns from Y for the objective
            %function
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            %Set up SDP
            % Y is the matrix representing all yields
            % J0 is the Choi matrix for the 0-photon component
            % J1 is the Choi matrix for the 1-photon component
            % In both Choi matrices Bob's system is first and Alice's
            % system is second to satisfy e.g. Y_0 = Tr[J_0 (F_y otimes (rho_x)^T)] 

            cvx_begin sdp quiet
                cvx_solver Mosek
                cvx_precision medium
                variable Y(m*(n_photon+1),n)
                variable J0(3,3) hermitian semidefinite
                variable J1(6,6) hermitian semidefinite
                maximize Objrow*Y*Objcolumn
                subject to
                    % 0 <= Y <=1
                    vec(Y) >= vec(zeros(m*(n_photon+1),n));
                    vec(Y) <= vec(ones(m*(n_photon+1),n));
                    
                    %0- and 1-photon component treated seperately with add.
                    %constraints
                    Y0 = M{1}*Y;
                    Y1 = M{2}*Y;
                    
                    %Partial trace of Choi matrix = id
                    norm(PartialTrace(J1,[1],[3,2])-eye(2)) <= choi_tolerance;
                    norm(trace(J0) - 1) <= choi_tolerance;
                    
                    %Usual decoy bounds rewritten as matrix times vector
                    for k = 1:n_decoy
                        pmu = Poisson(decoys(k),0:1:n_photon);
                        Pmu = kron(eye(m),pmu);
                        vec(Pmu*Y) >= vec(decoy_expecfull(:,:,k) - ones(m,n) * (1-sum(pmu)) - decoy_tolerance);
                        vec(Pmu*Y) <= vec(decoy_expecfull(:,:,k) + decoy_tolerance);                            
                    end

                    %Additional constraints for 0- and 1-photon components
                    %in terms of Choi matrices
                    for p = 1:m
                        for q = 1:n
                            norm(Y1(p,q) - trace(J1*kron(Fy{q},transpose(rhox{p})))) <= choi_tolerance;
                            norm(Y0(p,q) - trace(J0*kron(Fy{q},1))) <= choi_tolerance;
                        end
                    end
            cvx_end
            disp(cvx_status)
            Y1U(i+1,j) = Y(i*(n_photon+1)+2,j);
        end
    end

    %solve for lower bound
    disp("Solving for lower bounds")
    for i=0:m-1
        for j=1:n
            fprintf("\n Lower bounds solving SDP %d of %d ",n*i+(j-1)+1,m*n)

            
            %Select the correct rows and columns from Y for the objective
            %function
            Objcolumn = zeros(n,1);
            Objcolumn(j) = 1;
            
            Objrow = zeros(1,m*(n_photon+1));
            Objrow(i*(n_photon+1)+2) = 1;
            
            %Set up SDP
            % Y is the matrix representing all yields
            % J0 is the Choi matrix for the 0-photon component
            % J1 is the Choi matrix for the 1-photon component
            % In both Choi matrices Bob's system is first and Alice's
            % system is second to satisfy e.g. Y_0 = Tr[J_0 (F_y otimes (rho_x)^T)] 
            
            cvx_begin sdp quiet
                cvx_solver Mosek
                cvx_precision medium
                
                variable Y(m*(n_photon+1),n)
                variable J0(3,3) hermitian semidefinite
                variable J1(6,6) hermitian semidefinite
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
                    norm(PartialTrace(J1,[1],[3,2])-eye(2)) <= choi_tolerance;
                    norm(trace(J0) - 1) <= choi_tolerance;
                    
                    %Usual decoy bounds rewritten as matrix times vector
                    for k = 1:n_decoy
                        pmu = Poisson(decoys(k),0:1:n_photon);
                        Pmu = kron(eye(m),pmu);
                        vec(Pmu*Y) >= vec(decoy_expecfull(:,:,k) - ones(m,n) * (1-sum(pmu)) - decoy_tolerance);
                        vec(Pmu*Y) <= vec(decoy_expecfull(:,:,k) + decoy_tolerance);                            
                    end
                    
                    %Additional constraints for 0- and 1-photon components
                    %in terms of Choi matrices
                    for p = 1:m
                        for q = 1:n
                            norm(Y1(p,q) - trace(J1*kron(Fy{q},transpose(rhox{p})))) <= choi_tolerance;
                            norm(Y0(p,q) - trace(J0*kron(Fy{q},1))) <= choi_tolerance;
                        end
                    end
            cvx_end
            disp(cvx_status)
            Y1L(i+1,j) = Y(i*(n_photon+1)+2,j);
        end
    end
end