function [binsVar, valVar, diffVar] = standardize_hist(diffVar)
%% Computes normalized histogram of the data (error in angles in degrees)
%input diffVar: error in angles
%output binsVar: normalized count
%output valVar: the hist bins
%output diffVar: the updated error in angles
% Mona Zehni, 2019

diffVar(find(diffVar>180)) = diffVar(find(diffVar>180)) - 360;
diffVar(find(diffVar<-180)) = diffVar(find(diffVar<-180)) + 360;
[binsVar, valVar] = hist(diffVar, 200);

vals = linspace(-180, 180, 360);
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

