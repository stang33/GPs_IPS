function [output,T] = StandardPCGWithTEDITMATCHPAPER(multVectByA, b, multByPInv, errorTol, iters, t)

    %Run CG from paper setup for mPCG.

    x_0 = zeros(length(b),1);
    r_0 = b;% - multVectByA(u_0);
    z_0 = multByPInv(r_0);
    p_0 = z_0;



    x_jminus1 = x_0;
    r_jminus1 = r_0;
    p_jminus1 = p_0;
    z_jminus1 = z_0;

    alphas = zeros(t,1);
    betas = zeros(t,1);

    for j = 1 : iters

        v_j = multVectByA(p_jminus1);
        alpha_j = (r_jminus1' * z_jminus1) / (v_j' * p_jminus1);
        x_j = x_jminus1 + alpha_j * p_jminus1;
        r_j = r_jminus1 - alpha_j * v_j;
        z_j = multByPInv(r_j);

        %if norm(r_j, 2) < errorTol && j > t
        if j > t

            diagonal = zeros(t,1);

            diagonal(1) = 1 / alphas(1);
            for i = 2 : t
        
                diagonal(i) = 1 / alphas(i) + betas(i-1) / alphas(i-1);
        
            end
        
        
            offdiagonal = zeros(t,1);
            for i = 1 : t-1
        
                offdiagonal(i) = sqrt(betas(i)) / alphas(i);
        
            end
            
            %upperDiag = [zeros(1) offdiagonal(:t-1)];
            newMatrix = zeros(t+1,1);
            newMatrix(2:t+1,1) = offdiagonal;
            upperdiag = newMatrix(1:t);
            matrixOut = spdiags([offdiagonal diagonal upperdiag],-1:1,t,t);
            T = full(matrixOut);

            output = x_j;
            return
        end

        beta_j = (r_j' * z_j) / (r_jminus1' *  z_jminus1);
        p_j =  z_j + beta_j * p_jminus1;

        p_jminus1 = p_j;
        r_jminus1 = r_j;
        x_jminus1 = x_j;
        z_jminus1 = z_j;

        if j <= t
            alphas(j) = alpha_j;
            betas(j) = beta_j;
        end

    end



    diagonal = zeros(t,1);

    diagonal(1) = 1 / alphas(1);
    for i = 2 : t

        diagonal(i) = 1 / alphas(i) + betas(i-1) / alphas(i-1);

    end


    offdiagonal = zeros(t,1);
    for i = 1 : t-1

        offdiagonal(i) = sqrt(betas(i)) / alphas(i); %this was alphas(i)

    end
    
    %upperDiag = [zeros(1) offdiagonal(:t-1)];
    newMatrix = zeros(t+1,1);
    newMatrix(2:t+1,1) = offdiagonal;
    upperdiag = newMatrix(1:t);
    matrixOut = spdiags([offdiagonal diagonal upperdiag],-1:1,t,t);
    T = full(matrixOut);

    output = u_jminus1;


        