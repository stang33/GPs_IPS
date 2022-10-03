function order = EMM(ds, k, M)

    %k = num pts to order
    %ds = 1d dist vals
    %M boxes

    lt = ds(1);
    rt = ds(length(ds));

    %Make M boxes.
    boxLen = (rt - lt) / (M + 1);

    %boxes stores indices of the highest pt in the box
    boxes = zeros(M,1);

    split = lt + boxLen;
    curpt = 1;
    for i = 1 : M
    
        %Split into necessary box.
        while ds(curpt) < split
            curpt = curpt + 1;
        end

        boxes(i) = curpt - 1;

        split = split + boxLen;     

    end

    
    %Boxes now contains indices of splits.

    %Next, order boxes from 1 to M+1.

    orderInd = zeros(M,1);

    orderInd(1) = floor(M / 2);

    %Break into recursive alg.

    pt = EMMRecur(0, orderInd(1), false, true);

    orderedSplits = zeros(M,1);

    orderedSplits(1) = orderInd(1);
    %binPlace recall


    for i = 1 : depth

        %Go through list and use all splits.
        for i = 1 : M

            %

            if orderedSplits(i) ~= 0

                

            end

        end

    end






end