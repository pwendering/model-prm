function v = createRandomFluxVector(minFlux,maxFlux)
        l = numel(minFlux);
        v = (maxFlux-minFlux) .* rand(l,1) + minFlux;
end