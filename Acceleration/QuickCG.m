function output = QuickCG(x_0, multByA, b, errorTol, sigma, n, D, M, L, U_re, P_r, P_c, rhoVect, newRawIndices)

    %This method solves Ax = b through CG.
    %Note: multByA is a passed method that multiples by sparse A.
    %Possible further improvements via preconditioning.

    r_0 = b - multByA(x_0, sigma, n, D, M, L, U_re, P_r, P_c, rhoVect, newRawIndices);
    v_0 = r_0;

    x_k = x_0;
    r_k = r_0;
    v_k = v_0;

    while norm(r_k) > errorTol

        %5 step CG.
        %Single multiplication by A per iteration.
        tempAv_k = multByA(v_k, sigma, n, D, M, L, U_re, P_r, P_c, rhoVect, newRawIndices);
        t_k = dot(r_k, r_k) / dot(v_k, tempAv_k);
        x_kplus1 = x_k + t_k * v_k;
        r_kplus1 = r_k - t_k * tempAv_k;
        s_kplus1 = dot(r_kplus1, r_kplus1) / dot(r_k, r_k);
        v_kplus1 = r_kplus1 + s_kplus1 * v_k;

        %Update variables to store.
        x_k = x_kplus1;
        r_k = r_kplus1;
        v_k = v_kplus1;

    end

    output = x_k;


        