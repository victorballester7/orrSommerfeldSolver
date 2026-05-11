classdef eddyVisc

    methods(Static)

    
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
                            [~,D2] = math.f_chebdiff(Ny,2);
                            [~, D4] = math.cheb4c(Ny); %Clamped boundary condition on D4 
                            
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
                            
                          
                            [~,D2] = math.f_chebdiff(Ny,2);
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
            [~,D2] = math.f_chebdiff(Ny,2);
            [~, D4] = math.cheb4c(Ny); %Clamped boundary condition on D4 
            
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

        function [q,lambda] = modalbasis_s(alpha,beta,Ny,nmodes)

            [y,A]=math.f_chebdiff(Ny,2);
            [~,W]=math.clencurt(Ny-1);
            
            Dy = squeeze(A(:,:,1));
            Dy2 = squeeze(A(:,:,2));
            
            
            %% Stokes operator laplacian(u) - grad(p) = lambda*u
            II = eye(size(Dy)); ZZ = zeros(size(II));
            k = sqrt(alpha^2+beta^2);
            Lapl = Dy2-k^2*II;
            
            L = [Lapl ZZ ZZ -1i*alpha*II;
                ZZ Lapl ZZ -Dy;
                ZZ ZZ Lapl -1i*beta*II;
                1i*alpha*II Dy 1i*beta*II ZZ];
            F = [II ZZ ZZ ZZ;
                ZZ II ZZ ZZ;
                ZZ ZZ II ZZ;
                ZZ ZZ ZZ ZZ];
            
            % BCs
            L(1,:) = 0; L(1,1) = 1.; F(1,:) = 0; %u(pm1)=0
            L(Ny,:) = 0; L(Ny,Ny) = 1.; F(Ny,:) = 0;
            L(Ny+1,:) = 0; L(Ny+1,Ny+1) = 1.; F(Ny+1,:) = 0; %v(pm1)=0
            L(2*Ny,:) = 0; L(2*Ny,2*Ny) = 1.; F(2*Ny,:) = 0;
            L(2*Ny+1,:) = 0; L(2*Ny+1,2*Ny+1) = 1.; F(2*Ny+1,:) = 0; %v(pm1)=0
            L(3*Ny,:) = 0; L(3*Ny,3*Ny) = 1.; F(3*Ny,:) = 0;
            
            
            q = zeros(3*Ny,nmodes);
            
            [V,lambda]= eig(L,F); %Stokes modes
            lambda = diag(lambda);
            pos = find(isfinite(lambda));
            lambda = lambda(pos); V=V(:,pos);
            pos = find(abs(lambda)<1e6);
            lambda = lambda(pos); V=V(:,pos);
            
            [~,pos] = sort(abs(lambda),'ascend');
            lambda = lambda(pos);
            V = V(:,pos);
            
            lambda = lambda(1:nmodes);
            
            q = V(:,1:nmodes);
            if(max(lambda(1:nmodes)<0)>0) disp(['Evaluating Stokes modes, alpha,beta = ' num2str(alpha) num2str(beta)]); end
        end
         

        function [A,B,C] = NS_Operators_full(Ny,Re,kxh,kzh,uDat,nuDat)

            U = uDat(1:end,1);%2:end-1
            Up = uDat(1:end,2);
            % Upp = uDat(1:end,3);
            
            Nut = nuDat(1:end,1) - 1;%2:end-1
            NutP = nuDat(1:end,2);
            % NutPP = nuDat(1:end,3);
            
            k2 = kxh^2 + kzh^2;
            Iint = eye(Ny);
            Zint = zeros(Ny);
            [~,D] = math.f_chebdiff(Ny,2);
            % [~, D4] = cheb4c(Ny); 
            
            D1 = D(1:end,1:end,1); %2:end-1,2:end-1,1
            D2 = D(1:end,1:end,2);
            DW=D1;
            
            % clearvars D2
            Lap = D2 - k2*eye(Ny);
            
            i = sqrt(-1);
            AU=i*kxh*diag(U);
            
            L1 = [-(AU - Lap/Re-diag(Nut)*Lap/Re -diag(NutP)*D1/Re), -diag(Up)+i*kxh*diag(NutP)/Re, -Zint, -i*kxh*Iint ];
            L2 = [-Zint, -(AU - Lap/Re-diag(Nut)*Lap/Re -diag(NutP)*D1/Re)+diag(NutP)*D1/Re, -Zint , -D1 ];
            L3 = [-Zint, i*kzh*diag(NutP)/Re, -(AU - Lap/Re-diag(Nut)*Lap/Re -diag(NutP)*D1/Re), -i*kzh*Iint ];
            L4 = [i*kxh*Iint, D1, i*kzh*Iint, Zint];
            L  = [L1; L2; L3; L4];
            %%%% compute omega matrix
            AI = [Iint,Zint,Zint,Zint;Zint,Iint,Zint,Zint;Zint,Zint,Iint,Zint;Zint,Zint,Zint,Zint];
            
            %%%%% Derive Neumann boundary conditions
            %%%%% dp/dy(y=-1 & y=1) = 0 (Neumann BC) 
            b1 = [DW(1,1) , DW(1,Ny)];   bn = [DW(Ny,1) , DW(Ny,Ny)]; 
            d1 = DW(1, 2:Ny-1);          dn = DW(Ny, 2:Ny-1);
            Pb = [b1 ; bn];              Pd = [d1 ; dn];
            P  = -(Pb^-1)*Pd;
            
            %%%%% Impose boundary conditions (BC) Dirichelet
            %%%%% Remove outer row and column to eliminate singular point
            LLNS = L;
            LLNS(1:Ny:4*Ny,:) = [];
            LLNS(Ny-1:Ny-1:4*(Ny-1),:) = [];
            LLNS(:,1:Ny:4*Ny) = [];
            LLNS(:,Ny-1:Ny-1:4*(Ny-1)) = [];
            
            ALNS = AI;
            ALNS(1:Ny:4*Ny,:) = [];
            ALNS(Ny-1:Ny-1:4*(Ny-1),:) = [];
            ALNS(:,1:Ny:4*Ny) = [];
            ALNS(:,Ny-1:Ny-1:4*(Ny-1)) = [];
            
            %%%% Impose boundary conditions (BC) Neumann on pressure
            %%%% q' = [u', v' , w' , p']
            Pc1 =  [L(2:Ny-1,3*Ny+1),      L(2:Ny-1,4*Ny)];
            LLNS(1:Ny-2, 3*(Ny-2)+1:4*(Ny-2)) =  LLNS(1:Ny-2, 3*(Ny-2)+1:4*(Ny-2))+ Pc1*P; 
            
            Pc2 =  [L(2+Ny:2*Ny-1,3*Ny+1),  L(2+Ny:2*Ny-1,4*Ny)];
            LLNS(1+(Ny-2):2*(Ny-2),   3*(Ny-2)+1:4*(Ny-2)) =  LLNS(1+(Ny-2):2*(Ny-2), 3*(Ny-2)+1:4*(Ny-2))+ Pc2*P;
            
            Pc3 = [L(2+2*Ny:3*Ny-1,3*Ny+1), L(2+2*Ny:3*Ny-1,4*Ny)];
            LLNS(1+2*(Ny-2):3*(Ny-2), 3*(Ny-2)+1:4*(Ny-2)) =  LLNS(1+2*(Ny-2):3*(Ny-2), 3*(Ny-2)+1:4*(Ny-2))+ Pc3*P; 
            
            Pc4 = [L(2+3*Ny:4*Ny-1,3*Ny+1), L(2+3*Ny:4*Ny-1,4*Ny)];
            LLNS(1+3*(Ny-2):4*(Ny-2), 3*(Ny-2)+1:4*(Ny-2)) =  LLNS(1+3*(Ny-2):4*(Ny-2), 3*(Ny-2)+1:4*(Ny-2))+ Pc4*P;
            
            %%%%% compute envelope of maximum transient growth 
            [U,S]=eig(LLNS,ALNS);
            S=diag(S);   %compute eigenvalues
            ind = isfinite(S);
            S = S(ind);
            U=U(:,ind);
            %S = sort(S,'descend');
            
            %%%Plot calculated eigenvalue results
            LR = real(S);
            LI = imag(S);
            
            [~,order]=sort(LR,'descend');
            S=S(order);
            U=U(:,order);
            
            A=diag(S);
            C=U(1:3*(Ny-2),:);
            % C=U(:,:);
            B=pinv(C);
            end


        function q = modalbasis_eddyVisc(Ny,Re,kx,kz,uDat,nuDat,eddyFactor,nmodes)

                    [A,B,C,~,~] = eddyVisc.OS_SQ_Operators(Ny,Re,kx,kz,uDat,nuDat,eddyFactor);
                    [A_a,B_a,C_a,~] = eddyVisc.OS_SQ_adjointOperators(Ny,Re,kx,kz,uDat,nuDat,eddyFactor);
                    Yinf_nut = lyap(A,A_a,B*B_a);
                    Rff_nut = C * Yinf_nut* C_a;
                
                    [eigVec,eigVal] = eig(Rff_nut); 
                    eigVal = diag(eigVal);
                    [~,eigValSort] = sort(abs(eigVal),'descend');
                    % [~,eigValSort] = sort(-real(eigVal));
                
                    eigVec_Nut = eigVec(:,eigValSort);
                
                    % eigVec_Nut = C*eigVec_Nut;
                
                    Nint = Ny - 2;
                
                    q = eigVec_Nut;
                
                    qaux = zeros(3*Ny,nmodes);
                    qaux(2:Ny-1,:) = q(1:Nint,1:nmodes);
                    qaux(Ny+2:2*Ny-1,:) = q(Nint+1:2*Nint,1:nmodes);
                    qaux(2*Ny+2:3*Ny-1,:) = q(2*Nint+1:3*Nint,1:nmodes);
                
                    q = qaux;
                
                
                
        end

        function q = modalbasis_eddyVisc_c(Ny,Re,kx,kz,uDat,nuDat,eddyFactor,nmodes)

                    [A,B,C,~,~] = OS_SQ_Operators(Ny,Re,kx,kz,uDat,nuDat,eddyFactor);
                    [A_a,B_a,C_a,~] = OS_SQ_adjointOperators(Ny,Re,kx,kz,uDat,nuDat,eddyFactor);
                    Wc = lyap(A,A_a,B*B_a);
                
                    [eigVec,eigVal] = eig(Wc); 
                    eigVal = diag(eigVal);
                    [~,eigValSort] = sort(abs(eigVal),'descend');
                    % [~,eigValSort] = sort(-real(eigVal));
                
                    eigVec_Nut = eigVec(:,eigValSort);
                
                    eigVec_Nut = C*eigVec_Nut;
                
                    Nint = Ny - 2;
                
                    q = eigVec_Nut;
                
                    qaux = zeros(3*Ny,nmodes);
                    qaux(2:Ny-1,:) = q(1:Nint,1:nmodes);
                    qaux(Ny+2:2*Ny-1,:) = q(Nint+1:2*Nint,1:nmodes);
                    qaux(2*Ny+2:3*Ny-1,:) = q(2*Nint+1:3*Nint,1:nmodes);
                
                    q = qaux;
                
                
                
                end

        function [q,eigVal] = modalbasis_eddyVisc_eig(Ny,Re,kx,kz,uDat,nuDat,eddyFactor,nmodes)

                    [A,B,C,~,~] = OS_SQ_Operators(Ny,Re,kx,kz,uDat,nuDat,eddyFactor);
                    [A_a,B_a,C_a,~] = OS_SQ_adjointOperators(Ny,Re,kx,kz,uDat,nuDat,eddyFactor);
                    Yinf_nut = lyap(A,A_a,B*B_a);
                    Rff_nut = C * Yinf_nut* C_a;
                
                    [eigVec,eigVal] = eig(Rff_nut); 
                    eigVal = diag(eigVal);
                    [~,eigValSort] = sort(abs(eigVal),'descend');
                
                    eigVec_Nut = eigVec(:,eigValSort);
                                
                    Nint = Ny - 2;
                
                    q = eigVec_Nut;
                
                    qaux = zeros(3*Ny,nmodes);
                    qaux(2:Ny-1,:) = q(1:Nint,1:nmodes);
                    qaux(Ny+2:2*Ny-1,:) = q(Nint+1:2*Nint,1:nmodes);
                    qaux(2*Ny+2:3*Ny-1,:) = q(2*Nint+1:3*Nint,1:nmodes);
                
                    q = qaux;
                    eigVal = eigVal(1:nmodes);
                
                
        end


        function [uDat,nuDat, eta] = baseflowVaryParam_Couette_old(Ny,Re, U)

                    [eta, DM] = math.f_chebdiff(Ny, 2);

                    utau = eddyVisc.get_utau(U,Re);

                    profile = U/utau;

                    edVisc = Re./(DM(:,:,1)*profile);

                    nuDat(:,1) = edVisc;
                    nuDat(:,2) = DM(:,:,1)*edVisc;
                    nuDat(:,3) = DM(:,:,2)*edVisc;

                    uDat(:,1) = profile;
                    uDat(:,2) = DM(:,:,1)*profile;
                    uDat(:,3) = DM(:,:,2)*profile;

        end

        function [uDat_dns,nuDat_dns, eta] = baseflowVaryParam_Couette(Ny,Re, U)

            [eta,DM] = math.f_chebdiff(Ny,2);
            tau_w = DM(:,:,1)*U*(1/Re);
            tau_w = tau_w(1);
            utau = sqrt(tau_w);
            Re_tau = utau/(1/Re);

            dUdy_tau = U./utau;
            dUdy_tau = DM(:,:,1)*dUdy_tau;
            visc = Re_tau./dUdy_tau;
            
            nuDat_dns(:,1) = visc;
            nuDat_dns(:,2) = DM(:,:,1)*visc;
            nuDat_dns(:,3) = DM(:,:,2)*visc;
            
            uDat_dns(:,1) = U;
            uDat_dns(:,2) = DM(:,:,1)*uDat_dns(:,1);
            uDat_dns(:,3) = DM(:,:,2)*uDat_dns(:,1);

        end

        function [uDat,nuDat, eta] = baseflowVaryParam_Poi(Ny,Re,eddyFactor,aa,kk)
            [eta, DM] = math.f_chebdiff(Ny, 2);
            nutPL = 0.5*(1 + kk^2 * Re^2/9 * (1 - eta.^2).^2 .* (1 + 2*eta.^2).^2 .* (1 - exp((abs(eta)-1)*Re/aa)).^2)  .^0.5 - 0.5; %nut/nu
            nuTPL = 1 + nutPL;
            
            RHS = -Re * eta;
            MEANop = diag(nuTPL) * DM(:,:,1); 
            
            MEANop = MEANop(:,2:end-1);
            
            U = MEANop \ RHS;
            U = [0;U;0];
            
            DU = DM(:,:,1)*U;
            D2U = DM(:,:,2)*U;
            
            DnutPL = DM(:,:,1)*nutPL;
            D2nutPL = DM(:,:,2)*nutPL;
            
            uDat = [U,DU,D2U];
            nuDat = [1+eddyFactor*nutPL,DnutPL,D2nutPL];
        
        end

        function [uDat,nuDat, eta] = VerifybaseflowVaryParam_Couette(Ny,Re,eddyFactor,aa,kk)
                nu = 1;
                [eta, DM] = math.f_chebdiff(Ny, 2);
                nutPL = (nu/2)*(1 + kk^2 * Re^2/9 * (1 - eta.^2).^2 .* (1 + 2*eta.^2).^2 .* (1 - exp((abs(eta)-1)*Re/aa)).^2)  .^0.5 - nu/2; %nut/nu
                nuTPL = 1 + nutPL;
                
                
                arr = Re*ones(Ny,1);
                RHS = diag(nuTPL)\arr;
                
                DM1 = DM;
                % DM1(1,:,1) = 0;
                % DM1(1,1,1) = 1;
                % DM(end,:,1) = 0;
                % DM(end,end,1) = 1;
                DM(1,:) = 0;
                DM(1,1) = 1;
                
                RHS(1) = 14.3982;
                % RHS(end) = -1;

                U = DM(:,:,1)\RHS;
                % U = u;
                % U = u52;
                

                DU = DM1(:,:,1)*U;
                D2U = DM1(:,:,2)*U;
                
                DnutPL = DM1(:,:,1)*nutPL;
                D2nutPL = DM1(:,:,2)*nutPL;
                
                uDat = [U,DU,D2U];
                nuDat = [1+eddyFactor*nutPL,DnutPL,D2nutPL];


        end
        


        function utau = get_utau(u, Re)
            

            len = length(u);

            [~,DM] = math.f_chebdiff(len,2);
        
            dudy_w = DM(:,:,1)*u;
            dudy_w = dudy_w(1);
        
            tau_w = dudy_w/Re;
        
            utau = sqrt(tau_w);

        end

        function q = modalbasis_full(Ny,Re,kxh,kzh,uDat,nuDat)

            [A,B,C] = eddyVisc.NS_Operators_full(Ny,Re,kxh,kzh,uDat,nuDat);

            Rgg = lyap(A,A',B * B');
            Ruu = C * Rgg * C';
            Ruu = 0.5 * (Ruu + Ruu');
            
            [eigVec,eigVal] = eig(Ruu); 
            
            eigVal = diag(eigVal);
            ind = isfinite(eigVal);
            eigVal = eigVal(ind);
            [~,eigValSort] = sort(real(eigVal),'descend');%abs()
            eigVal = eigVal(eigValSort);
            eigVal = diag(eigVal);
            eigVec = eigVec(:,eigValSort);
            eigVec = real(eigVec);
            q = eigVec;

            Nint = Ny-2;

            quaux = [zeros(1,3*(Ny-2));q(1:Nint,:);zeros(1,3*(Ny-2))];
            qvaux = [zeros(1,3*(Ny-2));q(Nint+1:2*Nint,:);zeros(1,3*(Ny-2))];
            qwaux = [zeros(1,3*(Ny-2));q(2*Nint+1:end,:);zeros(1,3*(Ny-2))];

            q = [quaux;qvaux;qwaux];

         end

        function [A,B,C,D,M] = OS_SQ_Operators_w(Ny,Re,kx,kz,uDat,nuDat,eddyFactor)

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
            [~,D2] = math.f_chebdiff(Ny,2);
            [~, D4] = math.cheb4c(Ny); %Clamped boundary condition on D4 
            
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
            
          
            [~,D2] = math.f_chebdiff(Ny,2);
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
            % B = [-i*kx*(LapDir\D1Dir),-k2*(LapDir\Iint), -i*kz*(LapDir\D1Dir); i*kz*Iint, Zint, -i*kx*Iint];
            B11 = -i*kx*(LapDir\(diag(Nut)*D1Dir + diag(NutP)));
            B12 = -k2*(LapDir\Iint); 
            B13 = -i*kz*(LapDir\(diag(Nut)*D1Dir) + diag(NutP));
            B21 = i*kz*diag(Nut).*Iint;
            B22 = Zint;
            B23 = -i*kx*diag(Nut).*Iint;

            
            % B = [-i*kx*(LapDir\(D1Dir*Nut + NutP)),-k2*(LapDir\Iint),-i*kz*(LapDir\D1Dir*Nut + NutP);i*kz*Nut*Iint, Zint, -i*kx*Nut*Iint];
            B = [B11 B12 B13;B21 B22 B23];
            % A = [LapNeu\(Los+LosNut), Zint; Couple, Lsq + LsqNut]; 
            % B = [-i*kx*(LapNeu\D1Dir),-k2*(LapNeu\Iint), -i*kz*(LapNeu\D1Dir); i*kz*Iint, Zint, -i*kx*Iint];
            
            % u = Cq, defined from continuity and vorticity defn
            C = 1/k2 * [i*kx*D1Dir, -i*kz*Iint; k2*Iint, Zint; i*kz*D1Dir, i*kx*Iint];
            % q = Du, from vorticity defn only 
            D = [Zint, Iint, Zint; i*kz*Iint, Zint, -i*kx*Iint];
            
            
            M = 1/k2*[-LapDir, Zint; Zint, Iint];

        end

        function [A_a,B_a,C_a,D_a] = OS_SQ_adjointOperators_w(Ny,Re,kx,kz,uDat,nuDat,eddyFactor)
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
            [~,D2] = math.f_chebdiff(Ny,2);
            [~, D4] = math.cheb4c(Ny); %Clamped boundary condition on D4 
            
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
            B_a = 1/k2 * [i*kx*(diag(Nut)*D1Dir + diag(NutP)), -i*kz*Iint; k2*Iint, Zint; i*kz*(diag(Nut)*D1Dir + diag(NutP)), i*kx*Iint];

            C_a = [-i*kx*(LapDir\D1Dir),-k2*(LapDir\Iint), -i*kz*(LapDir\D1Dir); i*kz*Iint, Zint, -i*kx*Iint];
            D_a = 1/k2 * [Zint, -i*kz*Iint; -LapDir,Zint; Zint, i*kx*Iint]; 



            B11 = -i*kx*(LapDir\(diag(Nut)*D1Dir + diag(NutP)));
            B12 = -k2*(LapDir\Iint); 
            B13 = -i*kz*(LapDir\(diag(Nut)*D1Dir) + diag(NutP));
            B21 = i*kz*diag(Nut).*Iint;
            B22 = Zint;
            B23 = -i*kx*diag(Nut).*Iint;
           

            % C_a = [B11 B12 B13;
            %        B21 B22 B23];
        end


    end


    


end