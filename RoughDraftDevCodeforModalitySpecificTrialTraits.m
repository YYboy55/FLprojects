{gui, hdgs, mods, cohs, deltas }

mods = { 1 2 3}
cohs= {0, [.2,.7], .7};

numHdgGroups = num times heading set needs be repeated
I want this to be a setting that takes variable input, with each modality able to receive different num reps multiplier.

ves x {v}
vis x {vis}
comb x {comb}

         gui , ``hdgs, ````mods, `` cohs, ``` deltas,  multiplier
argin = {{0}, {-12 0 12}, {1 3}, {.2 .2 .7}, {0 3 -3}, {3, 1, 1}};
any(ismember(mods,1)) * multiplier{1} + 
any(ismember(mods,2)) * length(cohs{2}) * multiplier{2} + 
any(ismember(mods,3)) * length(cohs{3}) * length(deltas) * multiplier{3}



