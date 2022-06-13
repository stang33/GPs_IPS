% Output

fprintf(['Error in parameters are [%.2e' char(177) '%.2e] \n'],mean(max(abs(errorhyp))),std(max(abs(errorhyp))));


fprintf(['Error in phie are [%.2e' char(177) '%.2e], [%.2e' char(177) '%.2e] \n'],mean(errorphis(1,:)),std(errorphis(1,:)),mean(errorphis(2,:)),std(errorphis(2,:)));

fprintf(['Error in phia are [%.2e' char(177) '%.2e], [%.2e' char(177) '%.2e] \n'],mean(errorphis(3,:)),std(errorphis(3,:)),mean(errorphis(4,:)),std(errorphis(4,:)));


fprintf(['Error in training trajectories are [%.2e' char(177) '%.2e], [%.2e' char(177) '%.2e] \n'],mean(errortrajs_train(1,:)),std(errortrajs_train(1,:)),mean(errortrajs_train(3,:)),std(errortrajs_train(3,:)));

fprintf(['Error in testing trajectories are [%.2e' char(177) '%.2e], [%.2e' char(177) '%.2e] \n'],mean(errortrajs_test(1,:)),std(errortrajs_test(1,:)),mean(errortrajs_test(3,:)),std(errortrajs_test(3,:)));
