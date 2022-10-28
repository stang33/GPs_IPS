function tau_star = HutchinsonEst(MultByA, l, dimOfziVect)


    %Run STE.
    tau_star = 0;
    for i = 1 : l
    
        %Make a random vector of ones and zeros.
        z_i = 2 * randi(2,dimOfziVect,1) - 3 * ones(dimOfziVect,1);
        tau_star = tau_star + z_i' * MultByA(z_i);
    
    end
    
    tau_star = (1 / l) * tau_star;





end