function [spectrum_struct] = XrayTubeSpectrumTasmip(tube_voltage,varargin)
% function [spectrum,energies,spectrum_struct] = XrayTubeSpectrumTasmip(tube_voltage,varargin)
%    Compute x-ray tube photon flux spectrum. The spectrum is either the energy at the anode or
%      the photon number flux (in photons/millimeter^2) at a distance of 1 meter from the tube anode 
%       for an exposure of 1 milli-amp-second (mAS)
%    Each value spectrum(k) array is  the flux with energy k-0.5 to k+0.5 keV
%    inputs:
%      tube_voltage: in kilovolts (kV) must be less than or equal to 140
%      window_thickness: (optional = 0) The thickness in mm of the aluminum tube window
%      ripple: (optional = 0) The tube voltage ripple in percent (0 to 100)
%      number_spectrum: (optional = false). If present, the output is a number spectrum
%         else the default is the energy spectrum
%     clipzeros: (optional = false) if present clips off the leading zeros
%     and adjusts energies
%    outputs:
%      spectrum: (photons/mm^2@1 meter or Joules at the anode per energy step) 
%                   an array of spectrum values for photon energy 1:numel(spectrum) keV
%      energies: the x-ray photon energies corresponding to each entry in spectrum
%      spectrum_struct: a standard spectrum struct. See code for fields 
%   example of use
%     s = XrayTubeSpectrumTasmip(140,'window_thickness',3,'ripple',3.5);
%
%  Translated from Boone & Seibert's genspec1.c See:
%        An accurate method for computer-generating tungsten anode x-ray spectra from 30 to 140 kV.
%        Boone JM, Seibert JA.  Med Phys. 1997 Nov;24(11):1661-70.
%  REA 3/2/09

   % process the input arguments
nreqargs = 1;
assert(nargin>=nreqargs);
assert( (tube_voltage>=30) && (tube_voltage<=140));
almm = 0;
ripple = 0;
dynamic_range = 1e10; % ratio of max to min spectrum values allowed % (gil) was 1e4
number_spectrum = false;
clipzeros = false;
if(nargin>nreqargs)
  i=1;
  while(i<=size(varargin,2))
     switch lower(varargin{i})
     case 'number_spectrum';         number_spectrum = true;
     case 'clipzeros';         clipzeros = true;
     case 'window_thickness';      almm=varargin{i+1};    assert(isnumeric(almm)&& (almm>=0) );    i=i+1;
     case 'ripple';         ripple=varargin{i+1};    assert(isnumeric(ripple) && (ripple>=0) && (ripple<=100));  i=i+1;
    otherwise
      error('Unknown argument %s given',varargin{i});
     end
     i=i+1;
  end
end

dat = GetDataLoc;
if ripple ==0
  spectrum_number = dc_spectral_model(tube_voltage,almm,dat);
else
      % do the ripple calculation. Assume full wave rectify so compute over
      % sine from 0 to pi , actually (0+dtheta):pi
    nsteps = 20;
    theta = linspace(0,pi,nsteps+1);
    theta = theta(2:end); % do not use 0 angle
    ripple_voltage_amplitude = 0.01*ripple*tube_voltage;
       % kv_ripple is an array with the voltages at each phase angle step
    kv_ripple = tube_voltage - ripple_voltage_amplitude + ripple_voltage_amplitude*abs(sin(theta));
    spec_accum = zeros(150,1);
    for k = 1:nsteps
      spec = dc_spectral_model(kv_ripple(k),almm,dat);
      spec_accum(1:numel(spec)) =  spec_accum(1:numel(spec)) + spec;
    end
    nenergies_output = ceil(max(kv_ripple(:)));
    spectrum_number = spec_accum(1:nenergies_output)/nsteps;
end

if number_spectrum
  spectrum = spectrum_number;
else
  spectrum = 1.60217733e-16*4*pi*(1e3)^2*(spectrum_number.*(1:numel(spectrum_number))');
end
energies = (1:numel(spectrum))';
if clipzeros
        % clip off the 0 spectrum values
        % use dynamic_range
    minspec = max(spectrum)/dynamic_range;
    idx = find(abs(spectrum)>minspec,1,'first');
    spectrum = spectrum(idx:end);
    spectrum_number = spectrum_number(idx:end);
    energies = energies(idx:end);
end
if nargout == 1
    spectrum_struct.energies = energies(:);
    spectrum_struct.dE = 1;
    spectrum_struct.specnum = spectrum_number(:);
    spectrum_struct.Nphotons = trapz( spectrum_struct.energies ,  spectrum_struct.specnum);
    spectrum_struct.specnum_norm = spectrum_struct.specnum/spectrum_struct.Nphotons;
    spectrum_struct.spectrum = spectrum; % energy spectrum density
    spectrum_struct.specegy = 1.60217733e-16*4*pi*(1e3)^2*(spectrum_number.*(1:numel(spectrum_number))');
    spectrum_struct.Ebar = sum(spectrum_struct.energies.*spectrum_struct.specnum_norm);
    spectrum_struct.E2bar = sum(spectrum_struct.energies.^2.*spectrum_struct.specnum_norm);
end
  
% ------------ Local Functions ---------------------------

function spec = dc_spectral_model(kvolt,almm,dat)
% Same as Boone&Seibert function
% Pass it data in dat struct. See GetDataLoc function
n = round(kvolt)+2; % number of values for spectrum, 1 keV increments 
   % c code-- returns kvolt+3 values. Not sure why
epoly = repmat([1 kvolt kvolt^2 kvolt^3],n,1);
spec = sum(dat.aa(1:n,:).*epoly,2);
if almm>0
  spec = spec.*exp(-0.1*2.7*almm*dat.mual(1:n));
end  
spec = spec(3:end); % return values from 1:kvolt. 


function dat = GetDataLoc
% returns the tasmip polynomial coefficients and aluminum mu
%  as a struct with fields aa and mual
dat.aa=[ ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[-2.470985e+000,+7.522494e-002,+1.601297e-004,+0.000000e+000]; ...
[-5.468520e+001,+2.825971e+000,-3.702585e-002,+1.685450e-004]; ...
[-1.149660e+002,+7.181666e+000,-1.041506e-001,+5.246942e-004]; ...
[-2.023117e+001,+7.523026e+000,-9.725916e-002,+5.262351e-004]; ...
[+3.440159e+002,+3.179575e+000,-3.306927e-002,+3.115530e-004]; ...
[+5.493292e+002,+1.507932e+001,-9.656648e-002,+5.142380e-004]; ...
[+1.032546e+003,+2.793458e+001,-1.517779e-001,+6.491026e-004]; ...
[+1.056836e+003,+5.305293e+001,+4.006209e-002,-7.164506e-004]; ...
[+1.098845e+003,+8.295003e+001,+3.061647e-001,-2.617126e-003]; ...
[+4.957978e+002,+1.470037e+002,-1.102818e-001,-1.354507e-003]; ...
[-1.437833e+002,+2.229100e+002,-6.206306e-001,+1.896847e-004]; ...
[-1.106664e+003,+2.770497e+002,-5.743618e-001,-1.066210e-003]; ...
[-2.281766e+003,+3.424422e+002,-5.793318e-001,-2.303580e-003]; ...
[-5.591722e+003,+4.724134e+002,-1.429958e+000,+5.049076e-004]; ...
[-9.340535e+003,+6.186368e+002,-2.407872e+000,+3.701711e-003]; ...
[-1.406504e+004,+7.760495e+002,-3.430400e+000,+6.646413e-003]; ...
[-1.920322e+004,+9.418671e+002,-4.544806e+000,+9.920156e-003]; ...
[-2.515954e+004,+1.130912e+003,-5.997636e+000,+1.441550e-002]; ...
[-3.151928e+004,+1.331120e+003,-7.556880e+000,+1.925802e-002]; ...
[-3.165938e+004,+1.293120e+003,-6.625241e+000,+1.593667e-002]; ...
[-3.197696e+004,+1.259429e+003,-5.721722e+000,+1.269609e-002]; ...
[-3.150203e+004,+1.213018e+003,-4.995401e+000,+1.068630e-002]; ...
[-3.404540e+004,+1.273283e+003,-5.440755e+000,+1.275048e-002]; ...
[-3.525747e+004,+1.267165e+003,-5.052590e+000,+1.140252e-002]; ...
[-3.659796e+004,+1.264495e+003,-4.698218e+000,+1.017435e-002]; ...
[-3.935522e+004,+1.325721e+003,-5.260133e+000,+1.251165e-002]; ...
[-4.239447e+004,+1.396684e+003,-5.961586e+000,+1.539180e-002]; ...
[-4.505477e+004,+1.445302e+003,-6.324550e+000,+1.657817e-002]; ...
[-4.807436e+004,+1.506528e+003,-6.841015e+000,+1.832282e-002]; ...
[-4.772176e+004,+1.455009e+003,-6.183720e+000,+1.609850e-002]; ...
[-4.687265e+004,+1.383587e+003,-5.296423e+000,+1.305337e-002]; ...
[-4.534002e+004,+1.304458e+003,-4.458635e+000,+1.029127e-002]; ...
[-4.729671e+004,+1.337299e+003,-4.768113e+000,+1.129840e-002]; ...
[-4.592165e+004,+1.239852e+003,-3.651701e+000,+7.505117e-003]; ...
[-4.417617e+004,+1.131552e+003,-2.422704e+000,+3.340713e-003]; ...
[-4.975325e+004,+1.307914e+003,-4.490898e+000,+1.093279e-002]; ...
[-5.613191e+004,+1.511968e+003,-6.875300e+000,+1.962943e-002]; ...
[-5.524074e+004,+1.421870e+003,-5.669106e+000,+1.487642e-002]; ...
[-5.449938e+004,+1.337319e+003,-4.527925e+000,+1.035718e-002]; ...
[-5.884185e+004,+1.478833e+003,-6.293272e+000,+1.687622e-002]; ...
[-6.310984e+004,+1.616216e+003,-8.009326e+000,+2.321589e-002]; ...
[-5.995594e+004,+1.496680e+003,-6.906032e+000,+1.977848e-002]; ...
[-5.964100e+004,+1.456697e+003,-6.534316e+000,+1.853666e-002]; ...
[-6.132553e+004,+1.489142e+003,-6.956800e+000,+2.005068e-002]; ...
[-6.304895e+004,+1.522434e+003,-7.390895e+000,+2.161122e-002]; ...
[-5.994340e+004,+1.380871e+003,-5.839743e+000,+1.619943e-002]; ...
[-5.610868e+004,+1.218272e+003,-4.092096e+000,+1.018410e-002]; ...
[-1.825729e+004,-1.382119e+002,+9.557819e+000,-2.140051e-002]; ...
[+2.220017e+004,-1.568661e+003,+2.389806e+001,-5.505689e-002]; ...
[+5.501707e+004,-2.721157e+003,+3.527805e+001,-8.047399e-002]; ...
[+8.922944e+004,-3.915854e+003,+4.704985e+001,-1.070557e-001]; ...
[+2.104991e+004,-1.557364e+003,+2.321886e+001,-5.134972e-002]; ...
[-5.076517e+004,+9.032211e+002,-1.579828e+000,+7.306299e-003]; ...
[-6.030789e+004,+1.202068e+003,-4.552311e+000,+1.419530e-002]; ...
[-6.984994e+004,+1.499854e+003,-7.513087e+000,+2.103801e-002]; ...
[-7.108636e+004,+1.507313e+003,-7.472137e+000,+2.024801e-002]; ...
[-7.327537e+004,+1.540893e+003,-7.689933e+000,+2.028554e-002]; ...
[-3.161176e+004,+1.297773e+002,+6.392479e+000,-1.693738e-002]; ...
[+1.036295e+004,-1.288012e+003,+2.051981e+001,-5.423905e-002]; ...
[-4.132485e+004,+4.420904e+002,+2.448595e+000,+2.202247e-005]; ...
[-9.983141e+004,+2.351143e+003,-1.722188e+001,+5.896824e-002]; ...
[-8.345827e+004,+1.820261e+003,-1.140761e+001,+3.474510e-002]; ...
[-6.038053e+004,+1.099142e+003,-3.836391e+000,+5.215208e-003]; ...
[-7.332230e+004,+1.472738e+003,-7.481134e+000,+1.644730e-002]; ...
[-8.866886e+004,+1.911744e+003,-1.172736e+001,+2.948703e-002]; ...
[-8.906282e+004,+1.903695e+003,-1.166640e+001,+2.953372e-002]; ...
[-9.122084e+004,+1.949906e+003,-1.212404e+001,+3.119028e-002]; ...
[-9.195919e+004,+1.956641e+003,-1.222022e+001,+3.155684e-002]; ...
[-9.393503e+004,+1.997570e+003,-1.264453e+001,+3.294245e-002]; ...
[-9.460591e+004,+1.985575e+003,-1.240631e+001,+3.188458e-002]; ...
[-9.465909e+004,+1.947305e+003,-1.191912e+001,+3.005542e-002]; ...
[-1.054958e+005,+2.287738e+003,-1.546565e+001,+4.192772e-002]; ...
[-1.128820e+005,+2.523280e+003,-1.806383e+001,+5.099440e-002]; ...
[-5.652375e+004,+8.460812e+002,-1.890296e+000,+0.000000e+000]; ...
[-6.253113e+004,+9.546213e+002,-2.421458e+000,+0.000000e+000]; ...
[-6.063249e+004,+9.093265e+002,-2.222830e+000,+0.000000e+000]; ...
[-5.839087e+004,+8.581494e+002,-1.999379e+000,+0.000000e+000]; ...
[-6.177439e+004,+9.096954e+002,-2.219623e+000,+0.000000e+000]; ...
[-6.551339e+004,+9.674375e+002,-2.466158e+000,+0.000000e+000]; ...
[-6.482105e+004,+9.463755e+002,-2.384063e+000,+0.000000e+000]; ...
[-6.396586e+004,+9.225355e+002,-2.290526e+000,+0.000000e+000]; ...
[-5.976377e+004,+8.384694e+002,-1.918134e+000,+0.000000e+000]; ...
[-5.483239e+004,+7.418415e+002,-1.492676e+000,+0.000000e+000]; ...
[-5.545914e+004,+7.392220e+002,-1.466754e+000,+0.000000e+000]; ...
[-5.191874e+004,+6.677125e+002,-1.159438e+000,+0.000000e+000]; ...
[-5.337262e+004,+6.864440e+002,-1.248563e+000,+0.000000e+000]; ...
[-5.499713e+004,+7.080823e+002,-1.349865e+000,+0.000000e+000]; ...
[-6.109855e+004,+8.103042e+002,-1.805236e+000,+0.000000e+000]; ...
[-6.780313e+004,+9.224389e+002,-2.301017e+000,+0.000000e+000]; ...
[-6.463570e+004,+8.536160e+002,-1.980542e+000,+0.000000e+000]; ...
[-6.142322e+004,+7.841977e+002,-1.658250e+000,+0.000000e+000]; ...
[-6.542573e+004,+8.551263e+002,-1.999140e+000,+0.000000e+000]; ...
[-6.850218e+004,+9.104404e+002,-2.275249e+000,+0.000000e+000]; ...
[-6.775178e+004,+8.733046e+002,-2.050653e+000,+0.000000e+000]; ...
[-5.670986e+004,+6.717305e+002,-1.174642e+000,+0.000000e+000]; ...
[-6.431161e+004,+7.982173e+002,-1.730212e+000,+0.000000e+000]; ...
[-7.284777e+004,+9.397040e+002,-2.345359e+000,+0.000000e+000]; ...
[-7.296366e+004,+9.370416e+002,-2.349089e+000,+0.000000e+000]; ...
[-7.251969e+004,+9.256901e+002,-2.318580e+000,+0.000000e+000]; ...
[-7.373791e+004,+9.387560e+002,-2.371741e+000,+0.000000e+000]; ...
[-7.522138e+004,+9.557057e+002,-2.440560e+000,+0.000000e+000]; ...
[-6.645010e+004,+8.129935e+002,-1.892077e+000,+0.000000e+000]; ...
[-5.391723e+004,+6.111141e+002,-1.110798e+000,+0.000000e+000]; ...
[-6.950106e+004,+8.381854e+002,-1.943843e+000,+0.000000e+000]; ...
[-7.656837e+004,+9.340291e+002,-2.272803e+000,+0.000000e+000]; ...
[-7.169818e+004,+8.562692e+002,-1.994058e+000,+0.000000e+000]; ...
[-6.307650e+004,+7.199495e+002,-1.490337e+000,+0.000000e+000]; ...
[-6.896102e+004,+8.014658e+002,-1.785938e+000,+0.000000e+000]; ...
[-7.948799e+004,+9.545463e+002,-2.356450e+000,+0.000000e+000]; ...
[-8.038940e+004,+9.603943e+002,-2.368062e+000,+0.000000e+000]; ...
[-8.186549e+004,+9.744751e+002,-2.411129e+000,+0.000000e+000]; ...
[-8.127234e+004,+9.784392e+002,-2.501457e+000,+0.000000e+000]; ...
[-6.447853e+004,+7.327550e+002,-1.638994e+000,+0.000000e+000]; ...
[-3.806982e+004,+3.131658e+002,+0.000000e+000,+0.000000e+000]; ...
[-3.797812e+004,+3.101094e+002,+0.000000e+000,+0.000000e+000]; ...
[-4.023389e+004,+3.255209e+002,+0.000000e+000,+0.000000e+000]; ...
[-4.280943e+004,+3.432826e+002,+0.000000e+000,+0.000000e+000]; ...
[-4.114666e+004,+3.272756e+002,+0.000000e+000,+0.000000e+000]; ...
[-3.925966e+004,+3.096545e+002,+0.000000e+000,+0.000000e+000]; ...
[+3.191650e+002,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[-4.425804e+004,+3.425401e+002,+0.000000e+000,+0.000000e+000]; ...
[+8.115607e+001,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[-3.867988e+004,+2.969811e+002,+0.000000e+000,+0.000000e+000]; ...
[+1.306709e+003,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+1.153422e+003,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+9.817065e+002,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+8.099662e+002,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+6.688839e+002,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+5.277812e+002,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+3.498336e+002,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+1.718605e+002,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
[+0.000000e+000,+0.000000e+000,+0.000000e+000,+0.000000e+000]; ...
    ];

dat.mual = [ ...
1.179256e+003;2.262452e+003;7.830310e+002;3.580597e+002; ... % deleted 0 as first entry since matlab has 1 based arrays
1.924010e+002;1.148566e+002;7.365411e+001;5.038457e+001;3.555154e+001; ...
2.604821e+001;1.967844e+001;1.533015e+001;1.226424e+001;9.731421e+000; ...
7.980188e+000;6.575550e+000;5.499624e+000;4.647433e+000;4.034430e+000; ...
3.423026e+000;3.016158e+000;2.609285e+000;2.332402e+000;2.058015e+000; ...
1.830070e+000;1.664097e+000;1.498124e+000;1.350686e+000;1.240501e+000; ...
1.130512e+000;1.044317e+000;9.581216e-001;8.951530e-001;8.321844e-001; ...
7.692157e-001;7.232670e-001;6.773182e-001;6.376852e-001;6.025698e-001; ...
5.674542e-001;5.431759e-001;5.188975e-001;4.946192e-001;4.703408e-001; ...
4.460625e-001;4.304834e-001;4.149043e-001;3.993702e-001;3.838705e-001; ...
3.683709e-001;3.577159e-001;3.470610e-001;3.364060e-001;3.257511e-001; ...
3.150961e-001;3.076372e-001;3.001783e-001;2.927194e-001;2.852604e-001; ...
2.778015e-001;2.719912e-001;2.661808e-001;2.603705e-001;2.560495e-001; ...
2.517284e-001;2.474074e-001;2.430864e-001;2.387654e-001;2.344444e-001; ...
2.301956e-001;2.270014e-001;2.238071e-001;2.206129e-001;2.174186e-001; ...
2.144315e-001;2.118988e-001;2.093662e-001;2.068335e-001;2.043009e-001; ...
2.017682e-001;1.999037e-001;1.980392e-001;1.961747e-001;1.943102e-001; ...
1.924456e-001;1.905811e-001;1.887166e-001;1.868521e-001;1.849876e-001; ...
1.832079e-001;1.819418e-001;1.806756e-001;1.794095e-001;1.781434e-001; ...
1.768772e-001;1.756111e-001;1.743449e-001;1.730788e-001;1.718126e-001; ...
1.705465e-001;1.695537e-001;1.685609e-001;1.675681e-001;1.665754e-001; ...
1.655826e-001;1.645898e-001;1.635970e-001;1.626042e-001;1.616394e-001; ...
1.606749e-001;1.599088e-001;1.591427e-001;1.583766e-001;1.576106e-001; ...
1.568445e-001;1.560784e-001;1.553123e-001;1.546269e-001;1.539698e-001; ...
1.533127e-001;1.526556e-001;1.519985e-001;1.513414e-001;1.506843e-001; ...
1.500266e-001;1.494990e-001;1.489714e-001;1.484438e-001;1.479162e-001; ...
1.473885e-001;1.468609e-001;1.463333e-001;1.458057e-001;1.452781e-001; ...
1.447505e-001;1.442229e-001;1.436953e-001;1.431676e-001;1.426400e-001; ...
1.421124e-001;1.416803e-001;1.412483e-001;1.408162e-001;1.403841e-001; ...
1.399521e-001;1.395200e-001;1.390879e-001;1.386559e-001;1.382238e-001; ...
1.378063e-001;1.374185e-001;1.370307e-001;1.366429e-001;1.362550e-001; ...
1.358672e-001;1.354794e-001;1.350916e-001;1.347038e-001;1.343159e-001; ...
1.339281e-001;1.335927e-001;1.332572e-001;1.329218e-001;1.325863e-001; ...
1.322509e-001;1.319154e-001;1.315800e-001;1.312445e-001;1.309091e-001; ...
1.305736e-001;1.302595e-001;1.299482e-001;1.296369e-001;1.293256e-001; ...
1.290143e-001;1.287030e-001;1.283918e-001;1.280805e-001;1.277692e-001; ...
1.274579e-001;1.271811e-001;1.269044e-001;1.266276e-001;1.263508e-001; ...
1.260741e-001;1.257973e-001;1.255206e-001;1.252438e-001;1.249670e-001; ...
1.246903e-001;1.244135e-001;1.241367e-001;1.238600e-001;1.235832e-001; ...
1.233064e-001;1.230297e-001;1.227529e-001;1.224761e-001;1.222038e-001; ...
1.219318e-001];
