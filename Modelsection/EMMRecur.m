function pt = EMMRecur(start, en, trueL, trueR)

    %trueL, R signal which are incorporated into data.

    pt = -1;

    if trueL && trueR

        pt = (start - en) / 2 + start;

    elseif trueL

        pt = ceil(((start - en) - 1) / 3) + start;

    else

        pt = floor((2*(start - en) + 1) / 3) + start;

    end








end