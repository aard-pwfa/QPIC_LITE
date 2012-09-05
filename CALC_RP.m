function [output_struct] = CALC_RP(input_struct)

% import standard SI constants
eval(['run ' pwd '/SI_consts.m']);

% Calc plasma parameters
output_struct.omega_p  = sqrt(input_struct.plasma_density*1e6*SI_e^2/(SI_em*SI_eps0));  % plasma frequency  [rad/s]
output_struct.lambda_p = 2*pi*SI_c*1e6/output_struct.omega_p;                           % plasma wavelength [um]
output_struct.k_p      = output_struct.omega_p/SI_c;                                    % plasma wavenumber [1/m]
output_struct.SD       = 1e6/output_struct.k_p;                                         % plasma skin depth [um]

% Transfer used parameters
output_struct.plasma_density = input_struct.plasma_density;

% Calc gamma if not specified
if input_struct.gamma == 0
   if input_struct.energy == 0
      error('Must specify either gamma or energy');
   end
   output_struct.gamma  = input_struct.energy/(input_struct.mass*SI_eM/1e3);
   output_struct.energy = input_struct.energy;
end

% Calc energy if not specified
if input_struct.energy == 0
   if input_struct.gamma == 0
      error('Must specify either gamma or energy');
   end
   output_struct.energy = input_struct.gamma*input_struct.mass*SI_eM/1e3;
   output_struct.gamma  = input_struct.gamma;
end

% Betatron frequecy, wavelength
output_struct.omega_b = sqrt(output_struct.plasma_density*1e6*SI_e^2/...
    (2*output_struct.gamma*input_struct.mass*SI_em*SI_eps0));
output_struct.lambda_b = 2*pi*SI_c*1e6/output_struct.omega_b;

% Calc beam size if beam_match == 1
if( input_struct.beam_match )
  output_struct.sigma_x = 1e3*sqrt(input_struct.emit_x/output_struct.k_p*sqrt(2/output_struct.gamma)); % [um]
  output_struct.sigma_y = 1e3*sqrt(input_struct.emit_y/output_struct.k_p*sqrt(2/output_struct.gamma)); % [um]
  output_struct.emit_x  = input_struct.emit_x; % [mm*mrad]
  output_struct.emit_y  = input_struct.emit_x; % [mm*mrad]
end

% Calc emittance if emit_match == 1
if( input_struct.emit_match )
  output_struct.emit_x  = 1e-6*output_struct.k_p*(input_struct.sigma_x)^2*sqrt(output_struct.gamma/2); % [mm*mrad]
  output_struct.emit_y  = 1e-6*output_struct.k_p*(input_struct.sigma_y)^2*sqrt(output_struct.gamma/2); % [mm*mrad]
  output_struct.sigma_x = input_struct.sigma_x; % [um]
  output_struct.sigma_y = input_struct.sigma_y; % [um]
end

% Calc sigma_z if z_match == 1
if( input_struct.z_match )
   output_struct.sigma_z = 1e6*sqrt(2)/output_struct.k_p; % [um]
else
   output_struct.sigma_z = input_struct.sigma_z;
end

% Calc beam density
output_struct.beam_density = (10000^3)*input_struct.N_particles/...
    ((2*pi)^(3/2)*output_struct.sigma_x*output_struct.sigma_y*output_struct.sigma_z);

% Calc N_b/N_p
output_struct.ratio = output_struct.beam_density/output_struct.plasma_density;

% Calc k_p*sigma_x,y,z
output_struct.kpsx = output_struct.k_p * output_struct.sigma_x/1e6;
output_struct.kpsy = output_struct.k_p * output_struct.sigma_y/1e6;
output_struct.kpsz = output_struct.k_p * output_struct.sigma_z/1e6;

% Calc peak gaussian current
output_struct.I_peak = 1e3*input_struct.N_particles*SI_e*SI_c/(sqrt(2*pi)*output_struct.sigma_z); % [kA]
output_struct.N_particles   = input_struct.N_particles;

% Calc max bubble radius
output_struct.R_bubble = 1e6*(1/0.84)*2/output_struct.k_p*sqrt(output_struct.I_peak/17); % [um]

% Determine simulation box size
Box_X_exact = max(input_struct.N_bubbleradius*output_struct.R_bubble, input_struct.N_bunchradius*output_struct.sigma_x);  % [um]
Box_Y_exact = max(input_struct.N_bubbleradius*output_struct.R_bubble, input_struct.N_bunchradius*output_struct.sigma_y);  % [um]
Box_Z_exact = max(input_struct.N_bunchlengths*output_struct.sigma_z,  input_struct.N_wavelengths*output_struct.lambda_p); % [um]

% Round to nearet micron
output_struct.Box_X = ceil(Box_X_exact/10)*10; % [um]
output_struct.Box_Y = ceil(Box_Y_exact/10)*10; % [um]
output_struct.Box_Z = ceil(Box_Z_exact/10)*10; % [um]

% Determine grid spacing in terms of skin depth
skin_frac = 1/20;

% Determine number of cells, must be at least 2^6 per dim
output_struct.INDX = max(6, floor( log(output_struct.Box_X*1e-6 * output_struct.k_p / skin_frac) / log(2) ));
output_struct.INDY = max(6, floor( log(output_struct.Box_Y*1e-6 * output_struct.k_p / skin_frac) / log(2) ));
output_struct.INDZ = max(6, floor( log(output_struct.Box_Z*1e-6 * output_struct.k_p / skin_frac) / log(2) ));

% Increase granularity if desired
output_struct.INDX = output_struct.INDX + input_struct.x_grain;
output_struct.INDY = output_struct.INDY + input_struct.y_grain;
output_struct.INDZ = output_struct.INDZ + input_struct.z_grain;

% Determine cell size
output_struct.Cell_X = output_struct.Box_X/(2^(output_struct.INDX)); % [um]
output_struct.Cell_Y = output_struct.Box_Y/(2^(output_struct.INDY)); % [um]
output_struct.Cell_Z = output_struct.Box_Z/(2^(output_struct.INDZ)); % [um]

% Do not allow transverse cell size to be smaller than longitudinal cell size
while output_struct.Cell_Z > output_struct.Cell_X
   output_struct.INDZ = output_struct.INDZ + 1;
   output_struct.Cell_Z = output_struct.Box_Z/(2^(output_struct.INDZ));
end

% Cell Fraction
output_struct.Frac_X = output_struct.Cell_X/output_struct.SD;
output_struct.Frac_Z = output_struct.Cell_Z/output_struct.SD;

% Place beam on even/odd cell
output_struct.Box_X = output_struct.Box_X + input_struct.even_odd;
output_struct.Box_Y = output_struct.Box_Y + input_struct.even_odd;

% Determine beam position
output_struct.beam_x_pos = round(output_struct.Box_X*input_struct.x_frac);
output_struct.beam_y_pos = round(output_struct.Box_Y*input_struct.y_frac);
output_struct.beam_z_pos = round(output_struct.Box_Z*input_struct.z_frac);

% Determine slice time
output_struct.DT = round(sqrt(2*output_struct.gamma)/10 / 1.5); % 1.5 is deceleration factor
if( input_struct.BEAM_EV )
  output_struct.TEND       = floor(input_struct.plasma_s_prop / (SI_c / output_struct.omega_p))+0.1;
  output_struct.DT_OUTPUT  = 10; % output data every n'th timestep
  output_struct.DT_RESTART = output_struct.DT_OUTPUT*4;
else
  output_struct.TEND       = output_struct.DT+0.1;
  output_struct.DT_OUTPUT  = 1; % output data every timestep
  output_struct.DT_RESTART = 0;
end

% copy over remaining variables
output_struct.plasma_PREION = input_struct.plasma_PREION;
output_struct.plasma_Z      = input_struct.plasma_Z;
output_struct.charge        = input_struct.charge;
output_struct.mass          = input_struct.mass;
output_struct.store_QEB_3D  = input_struct.store_QEB_3D;


% Dump values?
if 0
   output_struct
end
