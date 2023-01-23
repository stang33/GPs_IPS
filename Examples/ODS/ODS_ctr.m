function f = ODS_ctr(x,i,s)

% input: P_i is the bias of stubborn agent i, R^d 
% (c) XXXX

d = length(x);
k = 10;
P_a = 1;
P_b = -1;
P_c = 0;

switch s
    case 1

        if i == 1
            f = -k*(x-P_a*ones(d,1));
        else
            f = 0*ones(d,1);
        end

    case 2

        if i == 1
            f = -k*(x-P_a*ones(d,1));
        elseif i == 2
            f = -k*(x-P_b*ones(d,1));
        else
            f = 0*ones(d,1);
        end     
        
     case 3

        if i == 1
            f = -k*(x-P_a*ones(d,1));
        elseif i == 2
            f = -k*(x-P_b*ones(d,1));
        elseif i == 3
            f = -k*(x-P_c*ones(d,1));
        else
            f = 0*ones(d,1);
        end     
        

        
end

end