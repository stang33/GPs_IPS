function output = RunCGOnMatrixInitGuesser(multByA, B, multByPInv, errorTol, iters)

    %B is n x t, where n is the size of A.
    output = 0 * B;

    initGuess = zeros(size(B,1),1);

    %For each column of B, run standard CG.
    %Use preguessing of previous answer.
    %Faster than parallelization on reasonable size data.
    for i = 1 : size(B,2)

        initGuess = StandardPCGInitGuesser(multByA, B(:,i), multByPInv, errorTol, iters, initGuess);
        output(:,i) = initGuess;

    end


        