function tau_star = HutchPlusPlus(MultMatrixByA, l, dimOfziVect)

    columnsInRand = floor(l / 3);

    %Make two random matrices of +-1s.
    S = 2 * randi(2,dimOfziVect,columnsInRand) - 3 * ones(dimOfziVect,columnsInRand);
    G = 2 * randi(2,dimOfziVect,columnsInRand) - 3 * ones(dimOfziVect,columnsInRand);

    %[Q,~] = qr(MultMatrixByA(S),"econ");
    [Q,~] = qr(MultMatrixByA(S));
    
    tau_star = trace(Q' * MultMatrixByA(Q)) + (3 / l) * trace(G' * (eye(dimOfziVect) - Q * Q') * MultMatrixByA((eye(dimOfziVect) - Q * Q') * G));

end