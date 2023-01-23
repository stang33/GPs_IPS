function out = TurnVectorToMatrixMult(MultVectByA, M, rows)

    columns = size(M,2);

    out = zeros(rows,columns);

    for i = 1 : columns

        out(:,i) = MultVectByA(M(:,i));

    end


end