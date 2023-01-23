function output = StandardPCG(multVectByA, b, multByPInv, errorTol, iters, initGuess)

    %Run CG from paper setup for mPCG.

    u_0 = zeros(length(b),1);
    if ~(initGuess == 0)
        u_0 = initGuess;
    end
    r_0 = b - multVectByA(u_0);

    %Avoid 0 solutions returning NaN.
    if norm(r_0) < errorTol
        output = u_0;
        return
    end

    z_0 = multByPInv(r_0);
    d_0 = z_0;

    d_jminus1 = d_0;
    r_jminus1 = r_0;
    z_jminus1 = z_0;
    u_jminus1 = u_0;

    for j = 1 : iters

        v_j = multVectByA(d_jminus1);
        alpha_j = (r_jminus1' * z_jminus1) / (d_jminus1' * v_j);
        u_j = u_jminus1 + alpha_j * d_jminus1;
        r_j = r_jminus1 - alpha_j * v_j;

        if norm(r_j, 2) < errorTol
            output = u_j;
            return
        end

        z_j = multByPInv(r_j);
        beta_j = (z_j' * r_j) / (z_jminus1' * r_jminus1);
        d_j = z_j + beta_j * d_jminus1;

        d_jminus1 = d_j;
        r_jminus1 = r_j;
        z_jminus1 = z_j;
        u_jminus1 = u_j;

    end

    output = u_jminus1;


        