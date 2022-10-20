function learnInfo = GetYm(learnInfo,hyp)


X = learnInfo.X;
Y = learnInfo.Y;
dN = learnInfo.d*learnInfo.N*learnInfo.order;
L = length(X)/dN;
d = learnInfo.d;
N = learnInfo.N;
%v = learnInfo.v;
name = learnInfo.name;

Nsub = N;
dN1 = d*N;
sub = 1:N;


if strcmp('ODS',name)
    Ym = zeros(dN1*L,1);
    logk = hyp(4); % parameter in control term
    Pa = hyp(5); % parameter in control term
    Pb = hyp(6); % parameter in control term
    Pc = hyp(7); % parameter in control term
    for i = 1:L
        Ym((dN1*(i-1)+1):(dN1*(i-1)+d)) = Y((dN1*(i-1)+1):(dN1*(i-1)+d)) + exp(logk)*(X((dN1*(i-1)+1):(dN1*(i-1)+d)) - Pa);
        Ym((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) = Y((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) + exp(logk)*(X((dN1*(i-1)+d+1):(dN1*(i-1)+2*d)) - Pb);
        Ym((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) = Y((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) + exp(logk)*(X((dN1*(i-1)+2*d+1):(dN1*(i-1)+3*d)) - Pc);
    end

    learnInfo.Ym = Ym;


elseif strcmp('CSF',name)
    Ym = zeros(dN1*L,1);
    loga = hyp(6); % parameter in noncollective force
    logb = hyp(7);
    for i = 1:L
        for j = 1:Nsub
            Ym((dN1*(i-1)+1+d*(j-1)):(dN1*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(sub(j)-1)):(dN1*(i-1)+d*sub(j))) - (exp(loga) * X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)))*(1- norm(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))),2)^exp(logb)));
        end
    end
    learnInfo.Ym = Ym;



elseif strcmp('FM',name)
    Ym = zeros(dN1*L,1);
    loga = hyp(6); % parameter in noncollective force
    logb = hyp(7);
    for i = 1:L
        for j = 1:Nsub
            Ym((dN1*(i-1)+1+d*(j-1)):(dN1*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(sub(j)-1)):(dN1*(i-1)+d*sub(j))) - (exp(loga) - exp(logb)*sum(X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j))).^2))*X((dN*(i-1)+1+dN1+d*(sub(j)-1)):(dN*(i-1)+dN1+d*sub(j)));
        end
    end
    learnInfo.Ym = Ym;
   

    
else
    Ym = zeros(dN1*L,1);
    for i = 1:L
        for j = 1:Nsub
            Ym((dN1*(i-1)+1+d*(j-1)):(dN1*(i-1)+d*j)) = Y((dN1*(i-1)+1+d*(sub(j)-1)):(dN1*(i-1)+d*sub(j)));
        end
    end

    learnInfo.Ym = Ym;
  
end


end



