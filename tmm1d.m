function DAT = tmm1d(DEV,SRC)
% INPUT ARGUMENTS
% ================
%% DEVELOPED BY: Ibrahim Issah.
%  School: University of Eastern Finland
%  Degree : Photonics
%
% .er1 relative permittivity in reflection region
% .ur1 relative permeability in reflection region
% .er2 relative permittivity in transmission region
% .ur2 relative permeability in transmission region
%
% .ER array containing permittivity of each layer
% .UR array containing permeability of each layer
% .L array containing thickness of each layer
%
% SRC Source Parameters
%
% .lam0 free space wavelength
%
% .theta elevation angle of incidence (radians)
% .phi azimuthal angle of incidence (radians)
%
% .pte amplitude of TE polarization
% .ptm amplitude of TM polarization
%
% OUTPUT ARGUMENTS
% ================
% DAT Output Data
%
% .R Reflectance
% .T Transmittance

er1 = DEV{1};
ur1 = DEV{2};
er2 = DEV{3};
ur2 = DEV{4};
ER = DEV{5};
UR = DEV{6};
L = DEV{7};

lam0 = SRC{1};
theta = SRC{2};
phi = SRC{3};
pte = SRC{4};
ptm = SRC{5};
%%
k0 = 2*pi/lam0;
kx = sqrt(ur1*er1)*sin(theta)*cos(phi);
ky = sqrt(ur1*er1)*sin(theta)*sin(phi);
Qh = [kx*ky 1+ky^2; -(1+kx^2) -kx*ky];
Vh = -1i*Qh;

SG11 = zeros(2,2);
SG12 = eye(2,2);
SG21 = eye(2,2);
SG22 = zeros(2,2);

for i = 1:length(L)
    Q = 1/UR(i) * [kx*ky  UR(i)*ER(i)-kx^2; ky^2-UR(i)*ER(i)  -kx*ky];
    kz = sqrt(UR(i)*ER(i) - kx^2 - ky^2);
    OMEGA = 1i*kz*eye(2,2);
    V = Q * OMEGA^(-1);
    X = expm(OMEGA*k0*L(i));
    A = eye(2,2) + V^(-1)*Vh;
    B = eye(2,2) - V^(-1)*Vh;
    D = A - X*B*A^(-1)*X*B;
    S11 = D^(-1)*(X*B*A^(-1)*X*A - B);
    S12 = D^(-1)*X*(A - B*A^(-1)*B);
    S21 = S12;
    S22 = S11;
    S = star({SG11 SG12 SG21 SG22},{S11 S12 S21 S22});
    SG11 = S{1};
    SG12 = S{2};
    SG21 = S{3};
    SG22 = S{4};
end

Q = 1/ur1 * [kx*ky  ur1*er1-kx^2; ky^2-ur1*er1  -kx*ky];
kz = sqrt(ur1*er1 - kx^2 - ky^2);
OMEGA = 1i*kz*eye(2,2);
Vref = Q * OMEGA^(-1);
Aref = eye(2,2) + Vh^(-1)*Vref;
Bref = eye(2,2) - Vh^(-1)*Vref;
SR11 = -Aref^(-1)*Bref;
SR12 = 2*Aref^(-1);
SR21 = 0.5*(Aref - Bref*Aref^(-1)*Bref);
SR22 = Bref*Aref^(-1);

Q = 1/ur2 * [kx*ky  ur2*er2-kx^2; ky^2-ur2*er2  -kx*ky];
kz = sqrt(ur2*er2 - kx^2 - ky^2);
OMEGA = 1i*kz*eye(2,2);
Vtrn = Q * OMEGA^(-1);
Atrn = eye(2,2) + Vh^(-1)*Vtrn;
Btrn = eye(2,2) - Vh^(-1)*Vtrn;
ST11 = Btrn*Atrn^(-1);
ST12 = 0.5*(Atrn - Btrn*Atrn^(-1)*Btrn);
ST21 = 2*Atrn^(-1);
ST22 = -Atrn^(-1)*Btrn;

S = star({SR11 SR12 SR21 SR22},{SG11 SG12 SG21 SG22});
SG11 = S{1};
SG12 = S{2};
SG21 = S{3};
SG22 = S{4};
S = star({SG11 SG12 SG21 SG22},{ST11 ST12 ST21 ST22});
SG11 = S{1};
SG12 = S{2};
SG21 = S{3};
SG22 = S{4};

n_inc = sqrt(ur1*er1);
k_inc = k0*n_inc*[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
n_surface_normal = [0;0;1];
if theta == 0
    a_te = [0;1;0];
else
    a_te = cross(n_surface_normal,k_inc) / normest(cross(n_surface_normal,k_inc));
end
a_tm = cross(a_te,k_inc) / normest(cross(a_te,k_inc));
p = pte*a_te + ptm*a_tm;

Esrc = [p(1,1);p(2,1)];
Eref = SG11*Esrc;
Etrn = SG21*Esrc;

Eztrn = - ( kx*Etrn(1,1) + ky*Etrn(2,1) ) /kz;
kzref = sqrt(ur1*er1)*cos(theta);
Ezref = - ( kx*Eref(1,1) + ky*Eref(2,1) ) / kzref;

R = ( normest([Eref;Ezref]) / normest(p) )^2;
T = ( normest([Etrn;Eztrn]) / normest(p) )^2 * real((ur1/ur2)*(kz/kzref));
%%
DAT = {R T};