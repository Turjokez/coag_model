addpath(genpath(pwd))
results = runtests("tests");
disp(table(results))
assert(all([results.Passed]), "Some tests failed. Stop and debug.");