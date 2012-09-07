% QuickPIC Matlab rpinput generation example script
% S. Gessner Sep 07, 2012

clear all;

% import standard SI constants
eval(['run ' pwd '/SI_consts.m']);

% specify template
rpinput_template_file = [pwd '/rpinputs/rpinput_template'];

% specify output
datestr = date;
[day, monthstr]  = strtok(datestr,'-');
[month, yearstr] = strtok(monthstr,'-');
[year, eonstr]   = strtok(yearstr,'-');

rpinput_dir = [pwd '/rpinputs/' year '/' month '/' day '/'];
if ~exist(rpinput_dir,'dir')
    mkdir(rpinput_dir);
end

param_dir = [pwd '/params/' year '/' month '/' day '/'];
if ~exist(param_dir,'dir')
    mkdir(param_dir);
end

rpinput_output_name = 'test1';
rpinput_output_file = [rpinput_dir 'rpinput_' rpinput_output_name];

write = 1;
% check to see if you want to overwrite file
if exist(rpinput_output_file,'file')
   reply = input(['File ' rpinput_output_file ' exists. \n Do you want to overwrite? y/n '], 's');
   if (strcmp(reply,'n'))
      disp('Ok. Forget it.');
      return;
   end
end


% INPUT TO RPINPUT

% plasma parameters
input_struct.plasma_density = 5e16;      % /cm^3
input_struct.plasma_PREION  = 1;           % 0 : non-ionized plasma 1: pre-ionized plasma
input_struct.plasma_Z       = 3;           % atomic number of plasma gas

% beam parameters
input_struct.charge         = +1.0;         % -1 for electron, +1 for positron
input_struct.mass           = SI_eM/SI_eM; % Particle mass in units of electron mass
input_struct.N_particles    = 2.0e10;        % Number of beam particles
input_struct.gamma          = 40000;           % relativistic factor gamma, if 0 energy specified below
input_struct.energy         = 0;          % beam mean energy [GeV], if 0 use gamma to calculate energy
input_struct.sigma_x        = 60;        % Gaussian sigma_x [um]
input_struct.sigma_y        = 60;        % Gaussian sigma_y [um]
input_struct.sigma_z        = 14;        % Gaussian sigma_z [um]
input_struct.emit_x         = 50.0;        % normalized X emittance [mm*mrad i.e. 1e-6]
input_struct.emit_y         = 50.0;         % normalized Y emittance [mm*mrad i.e. 1e-6]
input_struct.beam_match     = 0;           % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
input_struct.emit_match     = 0;           % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
input_struct.z_match        = 0;           % 1: override sigma_z with sqrt(2)/k_p, 0: do nothing

% simulation parameters
input_struct.BEAM_EV        = 0;           % 0 : calc wake only, 1 : propagate and evolve beam
input_struct.plasma_s_prop  = 1.0;         % propagation length of the beam [m]
input_struct.N_wavelengths  = 2.;           % set box length by number of plasma wavelengths
input_struct.N_bunchlengths = 0;          % set box length by bunch lengths, box length will be larger of two options
input_struct.N_bubbleradius = 0;        % set box width by number of bubble radii
input_struct.N_bunchradius  = 5.1;        % set box width by number of bunch radii, box width will be larger of two options
input_struct.x_frac         = 0.5;         % place bunch as fraction of box width, 0 at start of box, 1 at end
input_struct.y_frac         = 0.5;         % place bunch as fraction of box width, 0 at start of box, 1 at end
input_struct.z_frac         = 0.20;        % place bunch as fraction of box length, 0 at start of box, 1 at end
input_struct.even_odd       = 0;           % 0 : place bunch on even cell, 1 : place bunch on odd cell
input_struct.x_grain        = 0;           % increase granularity in the x dimension by 2^x_grain
input_struct.y_grain        = 0;           % increase granularity in the y dimension by 2^y_grain
input_struct.z_grain        = 0;           % increase granularity in the z dimension by 2^z_grain

% diagnostic parameters
input_struct.store_QEB_3D   = 0;           % store full 3D beam phase space?

% RPINPUT CALCULATOR
[output_struct] = CALC_RP(input_struct);

% RPINPUT WRITER
if write
    WRITE_RP(rpinput_template_file, rpinput_output_file, output_struct);
    save([param_dir 'param_' rpinput_output_name '.mat'], 'output_struct');
end
%exit;
