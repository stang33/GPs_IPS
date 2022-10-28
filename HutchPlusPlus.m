function tau_star = HutchPlusPlus(MultVectByA, l, dimOfziVect)

    %columnsInRand = floor(l / 3);
    columnsInRand = floor(l / 2);

    MultMatrixByA = @(M) TurnVectorToMatrixMult(MultVectByA, M, dimOfziVect);

    %Make two random matrices of +-1s.
    S = 2 * randi(2,dimOfziVect,columnsInRand) - 3 * ones(dimOfziVect,columnsInRand);
    G = 2 * randi(2,dimOfziVect,columnsInRand) - 3 * ones(dimOfziVect,columnsInRand);



    %[Q,~] = qr(MultMatrixByA(S),"econ");
    [Q,~] = qr(MultMatrixByA(S));

%     size(Q)
%     size(MultMatrixByA(Q))
%     size(Q' * MultMatrixByA(Q))
% 
%     t1 = trace(Q' * MultMatrixByA(Q))
%     t2 = trace(G' * (eye(dimOfziVect) - Q * Q') * MultMatrixByA((eye(dimOfziVect) - Q * Q') * G))
    
    tau_star = trace(Q' * MultMatrixByA(Q)) + (3 / l) * trace(G' * (eye(dimOfziVect) - Q * Q') * MultMatrixByA((eye(dimOfziVect) - Q * Q') * G));

end