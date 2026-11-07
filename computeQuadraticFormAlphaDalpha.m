%% -------------------------------------------------------------------------
function phi0 = computeQuadraticFormAlphaDalpha(Dcomp, Acomp)
% Compute phi0_e = alpha_e^T * D_e * alpha_e for each element
% Dcomp: 21 x nElem (lower-tri components of D)
% Acomp:  6 x nElem (alpha vector: [11,22,33,23,13,12])
% Returns: nElem x 1 vector
%
% Voigt lower-tri mapping (21 components):
% C1=D11, C2=D21, C3=D31, C4=D41, C5=D51, C6=D61
% C7=D22, C8=D32, C9=D42, C10=D52, C11=D62
% C12=D33, C13=D43, C14=D53, C15=D63
% C16=D44, C17=D54, C18=D64
% C19=D55, C20=D65
% C21=D66
%
% alpha^T * D * alpha = sum_i sum_j alpha_i * D_ij * alpha_j

a1 = Acomp(1,:)'; % 1,1 component
a2 = Acomp(2,:)'; % 2,2
a3 = Acomp(3,:)'; % 3,3
a4 = Acomp(4,:)'; % 2,3
a5 = Acomp(5,:)'; % 1,3
a6 = Acomp(6,:)'; % 1,2

phi0 = ...
    Dcomp(1,:)'.*a1.*a1 + Dcomp(7,:)'.*a2.*a2 + Dcomp(12,:)'.*a3.*a3 ...
  + Dcomp(16,:)'.*a4.*a4 + Dcomp(19,:)'.*a5.*a5 + Dcomp(21,:)'.*a6.*a6 ...
  + 2*Dcomp(2,:)'.*a1.*a2 + 2*Dcomp(3,:)'.*a1.*a3 + 2*Dcomp(4,:)'.*a1.*a4 ...
  + 2*Dcomp(5,:)'.*a1.*a5 + 2*Dcomp(6,:)'.*a1.*a6 + 2*Dcomp(8,:)'.*a2.*a3 ...
  + 2*Dcomp(9,:)'.*a2.*a4 + 2*Dcomp(10,:)'.*a2.*a5 + 2*Dcomp(11,:)'.*a2.*a6 ...
  + 2*Dcomp(13,:)'.*a3.*a4 + 2*Dcomp(14,:)'.*a3.*a5 + 2*Dcomp(15,:)'.*a3.*a6 ...
  + 2*Dcomp(17,:)'.*a4.*a5 + 2*Dcomp(18,:)'.*a4.*a6 + 2*Dcomp(20,:)'.*a5.*a6;

end