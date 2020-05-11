clear
clc

degrees = pi/180;

nPoints = 1000;
lam0 = linspace(1e-6,2e-6,nPoints); %free space wavelength
Rvector = zeros(1,nPoints);
Tvector = zeros(1,nPoints);

theta = 0 * degrees; %elevation angle
phi = 0 * degrees; %azimuthal angle
pte = 1; %amplitude of TE polarization
ptm = 0; %amplitude of TM polarization

ur1 = 1; %permeability in the reflection region
er1 = 1; %permittivity in the reflection region
ur2 = 1; %permeability in the transmission region
er2 = 1; %permittivity in the transmission region
% L=[150 250 150 250 150 250 150 250 150]*1e-9;
% ER=[3.5 1 3.5 1 3.5 1 3.5 1 3.5];
% UR=[1 1 1 1 1 1 1 1 1];
UR = [ 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00]; %array of permeabilities in each layer
ER = [ 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41 2.25 4.41]; %array of permittivities in each layer
L =  [ 250  180  250  180  250  180  250  180  250  180  250  180  250  180  250  180  250  180  250  180]*1e-9; %array of the thickness of each layer

DEV = {er1,ur1,er2,ur2,ER,UR,L};

count = 1;
for i = 1:length(lam0)
    SRC = {lam0(i),theta,phi,pte,ptm};
    DAT = tmm1d(DEV,SRC);
    Rvector(1,count) = DAT{1};
    Tvector(1,count) = DAT{2};
    count = count + 1;
end

CON = Rvector + Tvector; %% Conservation ==1 no loss > 1 gain <1 loss


%% PLOT SECTION 

figure('Color', 'white'); 

subplot1=subplot(2,1,1);hold on;box on;
plot(lam0*1e6,Rvector,'r', 'Linewidth', 2.5);
plot(lam0*1e6,Tvector,'b', 'Linewidth', 2.5);
plot(lam0*1e6,CON,'k--', 'Linewidth', 2.5);
ylim([0, 1]); 
xlabel('Wavelength (\mum)', 'FontSize',16); 
ylabel('Response','FontSize',16);
title('Wavelength response of a Bragg grating','FontSize',16); 
legend('Reflectance', 'Transmittance', 'Conservation'); 
legend1 = legend(subplot1,'show');
set(legend1, 'FontSize', 14,...
    'Position',[0.90625 0.873164218958611 0.07734375 0.0507343124165553]);

subplot2= subplot(2,1,2);hold on;box on;
plot(lam0*1e6,10*log10(Rvector),'r', 'Linewidth', 2.5);
plot(lam0*1e6,10*log10(Tvector),'b', 'Linewidth', 2.5);
plot(lam0*1e6,10*log10(CON),'k--', 'Linewidth', 2.5);
ylim([-30 0])
xlabel('Wavelength (\mum)', 'FontSize',16); 
ylabel('Response(dB)','FontSize',16); 
title('Wavelength response of a Bragg grating (log scale)','FontSize',16);
legend('Reflectance', 'Transmittance', 'Conservation'); 
legend2 = legend(subplot2,'show');
set(legend2,'Fontsize', 14,...
    'Position',[0.9046875 0.399198931909212 0.07734375 0.0507343124165554]);