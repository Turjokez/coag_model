function [fS,fM,fL,tot,used_mode] = export_size_classes(d_cm, flux, varargin)
% Export fractions using CM units
% Small  < 0.05 cm   (500 µm)
% Medium 0.05–0.2 cm (500–2000 µm)
% Large  > 0.2 cm    (>2000 µm)
% Optional: mode label ('volume' or 'image') for reporting

used_mode = 'volume';
if nargin >= 3 && ~isempty(varargin{1})
    used_mode = lower(string(varargin{1}));
end

d_cm  = d_cm(:);
flux  = flux(:);

good = isfinite(d_cm) & isfinite(flux);
d_cm = d_cm(good);
flux = max(flux(good),0);

T = sum(flux);

maskS = d_cm < 0.05;
maskM = d_cm >= 0.05 & d_cm <= 0.2;
maskL = d_cm > 0.2;

S = sum(flux(maskS));
M = sum(flux(maskM));
L = sum(flux(maskL));

if T > 0
    fS = S/T;
    fM = M/T;
    fL = L/T;
else
    fS = 0; fM = 0; fL = 0;
end

tot = struct('small',S,'medium',M,'large',L,'total',T);
end
