%% -------------------------------------------------------------------------
function eta_e = computeEtaPerElement(M_T, deltaTe)
% Compute eta_e = integral (DeltaT)^2 dOmega for each element
% deltaTe: 8 x nElem matrix of element temperature changes
% For each element: eta_e = DeltaT_e^T * M_T * DeltaT_e
nElem = size(deltaTe, 2);
eta_e = zeros(nElem, 1);
for e = 1:nElem
    dTe = deltaTe(:, e);
    eta_e(e) = dTe' * M_T * dTe;
end
end