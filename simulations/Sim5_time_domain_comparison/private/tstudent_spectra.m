function phi = tstudent_spectra(F, tstudent_a, tstudent_b, tstudent_dof, nu0)

% phi2=normpdf(f,mu,sigma);
if nargin==1
    tstudent_a = 600;
    tstudent_b = 9;
    nu0           = 10;
end
% phiXi=Axi./((1.0+(f./Bxi).^2).^3.2);
% phiAlpha=Aalpha./((1.0+((f-Ualpha)./Balpha).^2).^60);
phi = tstudent_a./((1.0+((F-nu0)./tstudent_b).^2).^tstudent_dof);
end



