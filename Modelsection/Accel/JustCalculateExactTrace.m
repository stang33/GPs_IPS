function tau_star = JustCalculateExactTrace(MultByA, dimOfziVect)


    %Run STE.
    tau_star = 0;
    parfor i = 1 : dimOfziVect
    
        %Make the unit vector.
        z_i = zeros(dimOfziVect,1);
        z_i(i,1) = 1;
        tau_star = tau_star + z_i' * MultByA(z_i);
    
    end


end