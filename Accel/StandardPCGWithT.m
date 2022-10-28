function [output,T] = StandardPCGWithT(multVectByA, b, multByPInv, errorTol, iters, t)

    %Run CG from paper setup for mPCG.
    %9/30: Changed for div by 0.

    u_0 = zeros(length(b),1);
    %r_0 = multVectByA(u_0) - b;
    r_0 = b - multVectByA(u_0);

    %Avoid 0 solutions returning NaN.
    %Note: fval doesn't like this. If
    %it's 0, send a signal to skip?
%     if norm(r_0) < errorTol
%         output = u_0;
%         return
%     end


    z_0 = multByPInv(r_0);
    d_0 = z_0;

    d_jminus1 = d_0;
    r_jminus1 = r_0;
    z_jminus1 = z_0;
    u_jminus1 = u_0;

    alphas = zeros(t);
    betas = zeros(t);

    for j = 1 : iters

        v_j = multVectByA(d_jminus1);
        alpha_j = (r_jminus1' * z_jminus1) / (d_jminus1' * v_j);
        u_j = u_jminus1 + alpha_j * d_jminus1;
        r_j = r_jminus1 - alpha_j * v_j;

        if norm(r_j, 2) < errorTol && j > t

            diagonal = zeros(t,1);

            diagonal(1) = 1 / alphas(1);
            for i = 2 : t
                
                if abs(alphas(i)) > 10^(-6)
                    diagonal(i) = 1 / alphas(i) + betas(i-1) / alphas(i-1);
                else
                    diagonal(i) = 0;
                end

            end
        
        
            offdiagonal = zeros(t,1);
            for i = 1 : t-1
        
                if abs(alphas(i)) > 10^(-6) && betas(i) > 0
                    offdiagonal(i) = sqrt(betas(i)) / alphas(i);
                else
                    offdiagonal(i) = 0;
                end

            end
            
            %upperDiag = [zeros(1) offdiagonal(:t-1)];
            newMatrix = zeros(t+1,1);
            newMatrix(2:t+1,1) = offdiagonal;
            upperdiag = newMatrix(1:t);
            matrixOut = spdiags([offdiagonal diagonal upperdiag],-1:1,t,t);
            T = full(matrixOut);
            
            output = u_j;
            return
        end

        z_j = multByPInv(r_j);
        %beta_j = (z_j' * z_j) / (z_jminus1' * z_jminus1);
        beta_j = (z_j' * r_j) / (z_jminus1' * r_jminus1);
        %d_j = z_j - beta_j * d_jminus1;
        d_j = z_j + beta_j * d_jminus1;

        d_jminus1 = d_j;
        r_jminus1 = r_j;
        z_jminus1 = z_j;
        u_jminus1 = u_j;

        if j <= t
            alphas(j) = alpha_j;
            betas(j) = beta_j;
        end

    end



    diagonal = zeros(t,1);

    diagonal(1) = 1 / alphas(1);
    for i = 2 : t

        if abs(alphas(i)) > 10^(-6)
            diagonal(i) = 1 / alphas(i) + betas(i-1) / alphas(i-1);
        else
            diagonal(i) = 0;
        end

    end


    offdiagonal = zeros(t,1);
    for i = 1 : t-1

        if abs(alphas(i)) > 10^(-6) && betas(i) > 0
            offdiagonal(i) = sqrt(betas(i)) / alphas(i);
        else
            offdiagonal(i) = 0;
        end

    end
    
    %upperDiag = [zeros(1) offdiagonal(:t-1)];
    newMatrix = zeros(t+1,1);
    newMatrix(2:t+1,1) = offdiagonal;
    upperdiag = newMatrix(1:t);
    matrixOut = spdiags([offdiagonal diagonal upperdiag],-1:1,t,t);
    T = full(matrixOut);

    output = u_jminus1;

        