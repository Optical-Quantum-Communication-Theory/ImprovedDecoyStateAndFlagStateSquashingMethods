%helper function that performs decoy state analysis
function [Y1L,Y1U] = SDP_decoy_BB84(decoys,decoy_expecfull)
       
    %cut-off for decoy
    n_photon=10;

    %number of intensities used
    n_decoy=size(decoys,2);

    %m = #number of rows
    %n = #number of colums
    m = size(decoy_expecfull,1);
    n = size(decoy_expecfull,2);
    
    %empty vectors for upper and lower bounds
    Y1L = zeros(m,n);
    Y1U = zeros(m,n);
    
    %Possonian distribution
    Poisson=@(mu,n) exp(-mu).*mu.^n./factorial(n);
    
    %Tolerances for Choi matrices and decoy constraints
    decoy_tolerance = 1e-12;
    choi_tolerance = 1e-12;
    
    %Signal states
    vec2dH = [1;0];
    vec2dV = [0;1];
    vec2dP = 1/sqrt(2)*(vec2dH + vec2dV);
    vec2dM = 1/sqrt(2)*(vec2dH - vec2dV);
    
    %Signal states as density matrices
    rhoH = vec2dH*vec2dH';
    rhoV = vec2dV*vec2dV';
    rhoP = vec2dP*vec2dP';
    rhoM = vec2dM*vec2dM';
    
    %combine into a list
    rhox = {rhoH,rhoV,rhoP,rhoM};
    
    %Bob's POVM elements
    vec3dH = [1;0;0];
    vec3dV = [0;1;0];
    vec3dP = 1/sqrt(2)*(vec3dH + vec3dV);
    vec3dM = 1/sqrt(2)*(vec3dH - vec3dV);
    vac = [0;0;1];
    
    FH = 1/2*vec3dH*vec3dH';
    FV = 1/2*vec3dV*vec3dV';
    FP = 1/2*vec3dP*vec3dP';
    FM = 1/2*vec3dM*vec3dM';
    Fvac = vac*vac';
    
    %combine into a list
    Fy = {FH,FV,FP,FM,Fvac};
    
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
            
            cvx_begin sdp quiet
            cvx_precision medium
            cvx_solver Mosek
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
                    
                    %Partial trace of Choi matrices = id
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

            cvx_begin sdp quiet
            cvx_precision medium
            cvx_solver Mosek
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

                    %Partial trace of Choi matrices = id
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