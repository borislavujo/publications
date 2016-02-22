Kl = load('Kl'); % load matrix
vinred = load('vi'); % load of list of states to keep
Kln = allngt(Kl,vinred); % eliminate all remaining states
save Kln Kln
