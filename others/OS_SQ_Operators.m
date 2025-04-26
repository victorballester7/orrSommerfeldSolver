function [A,B,C,D,M] = OS_SQ_Operators(Ny,Re,kx,kz,uDat,nuDat,eddyFactor)

                    Nint = Ny-2;
                    
                    U = uDat(2:end-1,1);
                    Up = uDat(2:end-1,2);
                    Upp = uDat(2:end-1,3);
                    
                    Nut = nuDat(2:end-1,1)-1;
                    NutP = nuDat(2:end-1,2);
                    NutPP = nuDat(2:end-1,3);


                    k2 = kx^2 + kz^2;
                    Iint = eye(Nint);
                    Zint = zeros(Nint);
                    [~,D2] = chebdif(Ny,2);
                    [~, D4] = cheb4c(Ny); %Clamped boundary condition on D4 
                    
                    D1Dir = D2(2:end-1,2:end-1,1); %Dirichlet Boundary conditions
                    D2Dir = D2(2:end-1,2:end-1,2);
                    clearvars D2
                    LapDir = D2Dir - k2*eye(Nint);
                    
                    Lap2 = D4 - 2*k2 * D2Dir + k2^2*Iint; %Clamped BCs on D4, else Dirichlet;
                    i = sqrt(-1);
                    
                    Los = -i*kx * (diag(U)*LapDir - diag (Upp)) + 1/Re * Lap2; % -i kx ( U Del - U'') + 1/Re (Del)^2 
                    Lsq = -i*kx * diag(U) + 1/Re * LapDir; % -i kx U + 1/Re Del;
                    
                     
                    LosNut = diag(Nut) * Lap2/Re + 2*diag(NutP) * (LapDir * D1Dir)/Re + diag(NutPP)*(D2Dir + k2*Iint)/Re;  
                    % LosNut = Nabla^2 N_y - d/dy (Nabla dot N), where N is nonlinear term  from nu_t
                    LsqNut = diag(Nut) * LapDir/Re + diag(NutP)*D1Dir/Re;
                    % LsqNut = d/dz N_x - d/dx N_z
                    Couple = -i*kz*diag(Up);
                    
                  
                    [~,D2] = chebdif(Ny,2);
                    LapNeu = D2(:,:,2) - k2*eye(Ny);
                    BC1 = D2(:,:,1);
                    BC1 = BC1([1,Ny],:);
                    G = -BC1(:,[1,Ny])\BC1(:,[2:Ny-1]);
                    
                    % LapNeu = LapNeu(2:Ny-1,2:Ny-1) + LapNeu(2:Ny-1,[1,Ny])*G;
                    
                    %% Construct state space
                    % d/dt  q =    [Lap, 0]^-1  [Los+LosNut,    0        ]q  +     B f 
                    %              [0  , I]     [-ikzU     , Lsq + LsqNut]        
                    % A = [LapDir\(Los+eddyFactor * LosNut), Zint; Couple, Lsq + eddyFactor * LsqNut]; 
                    % B = [-i*kx*(LapDir\D1Dir),-k2*(LapDir\Iint), -i*kz*(LapDir\D1Dir); i*kz*Iint, Zint, -i*kx*Iint];
                    
                    A = [LapDir\(Los+eddyFactor * LosNut), Zint; Couple, Lsq + eddyFactor * LsqNut]; 
                    B = [-i*kx*(LapDir\D1Dir),-k2*(LapDir\Iint), -i*kz*(LapDir\D1Dir); i*kz*Iint, Zint, -i*kx*Iint];
                    
                    % A = [LapNeu\(Los+LosNut), Zint; Couple, Lsq + LsqNut]; 
                    % B = [-i*kx*(LapNeu\D1Dir),-k2*(LapNeu\Iint), -i*kz*(LapNeu\D1Dir); i*kz*Iint, Zint, -i*kx*Iint];
                    
                    % u = Cq, defined from continuity and vorticity defn
                    C = 1/k2 * [i*kx*D1Dir, -i*kz*Iint; k2*Iint, Zint; i*kz*D1Dir, i*kx*Iint];
                    % q = Du, from vorticity defn only 
                    D = [Zint, Iint, Zint; i*kz*Iint, Zint, -i*kx*Iint];
                    
                    
                    M = 1/k2*[-LapDir, Zint; Zint, Iint];

end

function [A_a,B_a,C_a,D_a] = OS_SQ_adjointOperators(Ny,Re,kx,kz,uDat,nuDat,eddyFactor)
    Nint = Ny-2;
    
    U = uDat(2:end-1,1);
    Up = uDat(2:end-1,2);
    Upp = uDat(2:end-1,3);
    
    Nut = nuDat(2:end-1,1)-1; 
    NutP = nuDat(2:end-1,2);
    NutPP = nuDat(2:end-1,3);
    
    k2 = kx^2 + kz^2;
    Iint = eye(Nint);
    Zint = zeros(Nint);
    [~,D2] = chebdif(Ny,2);
    [~, D4] = cheb4c(Ny); %Clamped boundary condition on D4 
    
    D1Dir = D2(2:end-1,2:end-1,1); %Dirichlet Boundary conditions
    D2Dir = D2(2:end-1,2:end-1,2);
    clearvars D2
    LapDir = D2Dir - k2*eye(Nint);
    
    Lap2 = D4 - 2*k2 * D2Dir + k2^2*Iint; %Clamped BCs on D4, else Dirichlet;
    i = sqrt(-1);
    
    Los_a = i*kx * (diag(U)*LapDir + 2*diag (Up)*D1Dir) + 1/Re * Lap2;
    Lsq_a = i*kx * diag(U) + 1/Re * LapDir;
    
    LosNut_a = diag(Nut) * Lap2/Re + 2*diag(NutP) * (LapDir * D1Dir)/Re + diag(NutPP)*(D2Dir + k2*Iint)/Re;  
    LsqNut_a = diag(Nut) * LapDir/Re + diag(NutP)*D1Dir/Re;
    
    A_a = [LapDir\(Los_a + eddyFactor*LosNut_a), -i*kz*(LapDir\diag(Up)); Zint, Lsq_a + eddyFactor*LsqNut_a];
    B_a = 1/k2 * [i*kx*D1Dir, -i*kz*Iint; k2*Iint, Zint; i*kz*D1Dir, i*kx*Iint];
    C_a = [-i*kx*(LapDir\D1Dir),-k2*(LapDir\Iint), -i*kz*(LapDir\D1Dir); i*kz*Iint, Zint, -i*kx*Iint];
    D_a = 1/k2 * [Zint, -i*kz*Iint; -LapDir,Zint; Zint, i*kx*Iint]; 

end
