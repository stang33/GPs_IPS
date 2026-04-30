function trajData = trajUnifNoiseMultiplicative( trajData, sigma )
% Adds multiplicative gaussian noise of relative amplitude sigma to trajectories
% Authored by Jinchao Feng and Sui Tang

%trajData = trajData+sigma*(randn(size(trajData)));
trajData = trajData + sigma .* abs(trajData) .* randn(size(trajData));


return