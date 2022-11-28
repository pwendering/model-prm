function v = createRandomFluxVector(minFlux, maxFlux)
%% v = createRandomFluxVector(minFlux, maxFlux)
% Creates a randomvector with between minimal and maximal concentrations
% for each element in minFlux/maxFlux.
% INPUT
%   double minFlux:         array of minimum metabolite concentrations
%   double maxFlux:         array of maximum metabolite concentrations
% OUTPUT
%   double v:               array of random concentrations within the
%                           limits of minFlux and maxFlux for every element
l = numel(minFlux);
v = (maxFlux-minFlux) .* rand(l,1) + minFlux;
end