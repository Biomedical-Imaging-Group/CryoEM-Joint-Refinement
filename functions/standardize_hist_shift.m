function [binsVar, valVar] = standardize_hist_shift(diffVar, range, numBin)
%% Computes normalized histogram of the data (error in shifts)
%input diffVar: error in shifts
%input range: range of the histogram support
%input numBin: number of bins
%output binsVar: normalized count
%output valVar: the hist bins
% Mona Zehni, 2019

vals = linspace(-range, range, numBin);
binsVar = zeros(size(vals));
for i=1:length(vals)-1
    if i==length(vals)
        binsVar(i) = length(find(diffVar>=vals(i)));
    end
    binsVar(i) = length(find(diffVar<vals(i+1) & diffVar>=vals(i)));
end
binsVar = binsVar/sum(binsVar);
valVar = vals;

end