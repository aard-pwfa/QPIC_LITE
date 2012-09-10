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

rpinput_output_name = 'test2';
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

% simulation parameters
input_struct.sim.BEAM_EV       = 1;           % 0 : calc wake only, 1 : propagate and evolve beam
input_struct.sim.prop          = 0.000768;         % propagation length of the beam [m]
input_struct.sim.DT            = 16.0;        % Delta T between beam pushes [1/omega_p]. If 0: use calc from formula
input_struct.sim.dump_freq     = 1;           % Dump frequency

% plasma parameters
input_struct.plasma.density    = 5e16;        % /cm^3
input_struct.plasma.charge     = -1.0;        % -1 for electron, +1 for positron
input_struct.plasma.mass       = SI_eM/SI_eM; % Particle mass in units of electron mass
input_struct.plasma.PREION     = 1;           % 0 : non-ionized plasma 1: pre-ionized plasma
input_struct.plasma.Z          = 3;           % atomic number of plasma gas

% beam parameters
input_struct.beam.charge       = +1.0;        % -1 for electron, +1 for positron
input_struct.beam.mass         = SI_eM/SI_eM; % Particle mass in units of electron mass
input_struct.beam.N_particles  = 2.0e10;      % Number of beam particles
input_struct.beam.gamma        = 40000;       % relativistic factor gamma, if 0 energy specified below
input_struct.beam.energy       = 0;           % beam mean energy [GeV], if 0 use gamma to calculate energy
input_struct.beam.sigma_x      = 60;          % Gaussian sigma_x [um]
input_struct.beam.sigma_y      = 60;          % Gaussian sigma_y [um]
input_struct.beam.sigma_z      = 14;          % Gaussian sigma_z [um]
input_struct.beam.emit_x       = 50.0;        % normalized X emittance [mm*mrad i.e. 1e-6]
input_struct.beam.emit_y       = 50.0;        % normalized Y emittance [mm*mrad i.e. 1e-6]
input_struct.beam.beam_match   = 0;           % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
input_struct.beam.emit_match   = 0;           % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
input_struct.beam.z_match      = 0;           % 1: override sigma_z with sqrt(2)/k_p, 0: do nothing

% size parameters
input_struct.size.Z_waves      = 2.;          % set box length by number of plasma wavelengths
input_struct.size.Z_bunches    = 0;           % set box length by bunch lengths
input_struct.size.X_bubbles    = 0;           % set box width by number of bubble radii
input_struct.size.X_bunches    = 5.1;         % set box width by number of bunch radii
input_struct.size.X_center     = 0.5;         % place bunch as fraction of box width, 0 at start of box, 1 at end
input_struct.size.Y_center     = 0.5;         % place bunch as fraction of box width, 0 at start of box, 1 at end
input_struct.size.Z_center     = 0.20;        % place bunch as fraction of box length, 0 at start of box, 1 at end
input_struct.size.x_grain      = 0;           % increase granularity in the x dimension by 2^x_grain
input_struct.size.y_grain      = 0;           % increase granularity in the y dimension by 2^y_grain
input_struct.size.z_grain      = 0;           % increase granularity in the z dimension by 2^z_grain

% diagnostic parameters
input_struct.diag.store_QEB_3D = 0;           % store full 3D beam phase space?

% RPINPUT CALCULATOR
param_struct = CALC_RP(input_struct);

% RPINPUT WRITER
if write
    WRITE_RP(rpinput_template_file, rpinput_output_file, param_struct);
    save([param_dir 'param_' rpinput_output_name '.mat'], 'param_struct');
end
%exit;
