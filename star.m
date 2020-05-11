function S = star(SA,SB)
% STAR Redheffer Star Product
%
% S = star(SA,SB)
%
% INPUT ARGUMENTS
% ================
% SA First Scattering Matrix
% .S11 S11 scattering parameter
% .S12 S12 scattering parameter
% .S21 S21 scattering parameter
% .S22 S22 scattering parameter
%
% SB Second Scattering Matrix
% .S11 S11 scattering parameter
% .S12 S12 scattering parameter
% .S21 S21 scattering parameter
% .S22 S22 scattering parameter
%
% OUTPUT ARGUMENTS
% ================
% S Combined Scattering Matrix
% .S11 S11 scattering parameter
% .S12 S12 scattering parameter
% .S21 S21 scattering parameter
% .S22 S22 scattering parameter

SG11 = SA{1};
SG12 = SA{2};
SG21 = SA{3};
SG22 = SA{4};
S11 = SB{1};
S12 = SB{2};
S21 = SB{3};
S22 = SB{4};
D = SG12*(eye(2,2) - S11*SG22)^(-1);
F = S21*(eye(2,2) - SG22*S11)^(-1);
SG11 = SG11 + D*S11*SG21;
SG12 = D*S12;
SG21 = F*SG21;
SG22 = S22 + F*SG22*S12;
S = {SG11 SG12 SG21 SG22};