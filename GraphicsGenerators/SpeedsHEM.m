
clear('all');
clc;

rho = GeometricSpace(0.005,997,1.01,1000);
T = linspace(290,650,200);
[rho,T] = meshgrid(rho,T);

P = Pressure(rho,T);

dPdrho = Pressure_DensityN(rho,T);

dPdT = Pressure_Temperature(rho,T);
dTdi = 1./InternalEnergy_TemperatureN(rho,T);
dPdi = 1./rho .* dPdT .* dTdi;

Discriminant = @(u) 4 * P .* dPdi + (u .* dPdi).^2 + 4 * rho.^2 .* dPdrho;
Speed        = @(u) (1+dPdi./rho/2)*u + sqrt(Discriminant(u))./(2*rho);
u = linspace(0,50,500);
speedMin = u;
speedMax = u;
speedAvg = u;
DiscrimMin = u;

for k = 1:length(u)
	DiscrimMin(k) = min(min(Discriminant(u(k))));
    speedMin(k)   = min(min(Speed(u(k))));
    speedMax(k)   = max(max(Speed(u(k))));
    speedAvg(k)   = mean(mean(Speed(u(k))));
end

semilogy(u,speedMin,u,speedMax,u,speedAvg,'LineWidth',2);
FontSize = 10;
FontName = 'CMU Serif';
set(gca,'FontSize',FontSize,'FontName',FontName);
legend('Minimum speed','Maximum speed','Average speed',4);
legend(gca,'boxoff');
xlabel('Material Speed [m/s]','FontSize',FontSize,'FontName',FontName);
ylabel('Acoustic Speed [m/s]','FontSize',FontSize,'FontName',FontName);
grid('on');
box('on');
MakePNG('AcousticSpeedVsMaterialSpeed');
CropOuterWhite(WhereWeAre(),'AcousticSpeedVsMaterialSpeed');



