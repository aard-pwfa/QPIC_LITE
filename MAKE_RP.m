% QuickPIC Matlab rpinput generation example script
% S. Gessner Sep 07, 2012

% clear all;

input_struct.sim_name = 'testing2';

% INPUT TO RPINPUT

% import standard SI constants
SI_consts;

% simulation parameters
input_struct.sim.BEAM_EV       = 0;           % 0 : calc wake only, 1 : propagate and evolve beam
input_struct.sim.prop          = 1;           % propagation length of the beam [m]
input_struct.sim.DT            = 0;           % Delta T between beam pushes [1/omega_p]. If 0: use calc from formula
input_struct.sim.dump_freq     = 10;          % Dump frequency
input_struct.sim.Use_Destroyer = 'false';     % indicate 'true' or 'false' here, to enable (disable) the particle destroyer

% plasma parameters
input_struct.plasma.density    = 1e17;        % /cm^3
input_struct.plasma.charge     = -1.0;        % -1 for electron, +1 for positron
input_struct.plasma.mass       = SI_eM/SI_eM; % Particle mass in units of electron mass
input_struct.plasma.PREION     = 1;           % 0 : non-ionized plasma 1: pre-ionized plasma
input_struct.plasma.Z          = 1;           % atomic number of plasma gas
input_struct.plasma.profile    = 1;           % 0: uniform plasma, 1: hollow channel plasma
input_struct.plasma.n_point    = 5;           % number of points used to create plasma profile
input_struct.plasma.radius     = 33.61;       % channel radius [um]
input_struct.plasma.width      = 3;           % annulus width [um]

% beam parameters
input_struct.beam.charge       = -1.0;        % -1 for electron, +1 for positron
input_struct.beam.mass         = SI_eM/SI_eM; % Particle mass in units of electron mass
input_struct.beam.N_particles  = 1.0e10;      % Number of beam particles
input_struct.beam.gamma        = 39139;       % relativistic factor gamma, if 0 energy specified below
input_struct.beam.energy       = 0;           % beam mean energy [GeV], if 0 use gamma to calculate energy
input_struct.Init_Routine      = 1;           % 1 for gaussian beam initialization, 5 for Twiss parameter beam initialization
input_struct.beam.sigma_x      = 10;          % Gaussian sigma_x [um]
input_struct.beam.sigma_y      = 10;          % Gaussian sigma_y [um]
input_struct.beam.s_waist      = 0.35;        % Waist position relative to start of the simulation [m]
                                              % If 0, then alpha and beta indicated below are used. 
                                              % If different from 0, alpha and beta are calculated from
                                              % sigma_x, sigma_y and s_waist.
input_struct.beam.alpha_x      = 0;           % Twiss parameter alpha_x
input_struct.beam.beta_x       = 0.5;         % Twiss parameter beta_x [m]
input_struct.beam.alpha_y      = 0;           % Twiss parameter alpha_y
input_struct.beam.beta_y       = 5;           % Twiss parameter beta_y [m]
input_struct.beam.emit_x       = 50;          % normalized X emittance [mm*mrad i.e. 1e-6]
input_struct.beam.emit_y       = 5;           % normalized Y emittance [mm*mrad i.e. 1e-6]
input_struct.beam.sigma_z      = 20;          % Gaussian sigma_z [um]
input_struct.beam.beam_match   = 0;           % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
input_struct.beam.emit_match   = 0;           % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
input_struct.beam.z_match      = 0;           % 1: override sigma_z with sqrt(2)/k_p, 0: do nothing

% size parameters
input_struct.size.Z_waves      = 3.3;         % set box length by number of plasma wavelengths
input_struct.size.Z_bunches    = 0;           % set box length by bunch lengths
input_struct.size.X_bubbles    = 15;          % set box width by number of bubble radii
input_struct.size.X_bunches    = 0.0;         % set box width by number of bunch radii
input_struct.size.X_center     = 0.5;         % place bunch as fraction of box width, 0 at start of box, 1 at end
input_struct.size.Y_center     = 0.5;         % place bunch as fraction of box width, 0 at start of box, 1 at end
input_struct.size.Z_center     = 0.25;        % place bunch as fraction of box length, 0 at start of box, 1 at end
input_struct.size.x_grain      = 0;           % increase granularity in the x dimension by 2^x_grain
input_struct.size.y_grain      = 0;           % increase granularity in the y dimension by 2^y_grain
input_struct.size.z_grain      = 0;           % increase granularity in the z dimension by 2^z_grain

% diagnostic parameters
input_struct.diag.store_QEB_3D = 0;           % store full 3D beam phase space?

% computational parameters
input_struct.comp.qpic_executable = 'qpic.e.twiss.0907'; % Select the qpic executable to be used
if input_struct.sim.BEAM_EV == 0
    input_struct.comp.evolution  = '.false.';
    input_struct.comp.num_stages = 1;
    input_struct.comp.mem        = 1024;
    input_struct.comp.tasks      = 8;
    input_struct.comp.run_time   = 1;
elseif input_struct.sim.BEAM_EV == 1
    input_struct.comp.evolution  = '.true.';
    input_struct.comp.num_stages = 1;
    input_struct.comp.mem        = 2048;
    input_struct.comp.tasks      = 64;
    input_struct.comp.run_time   = 10;        % Amount of computer time to run sim for in hours
end




% CREATE DIRECTORIES AND VARIOUS VARIABLES

% specify template
rpinput_template_file = [pwd '/rpinputs/rpinput_template'];

% specify output
date_dir = GET_DATE_DIR;

rpinput_dir = [pwd '/rpinputs/' date_dir];
if ~exist(rpinput_dir,'dir')
    mkdir(rpinput_dir);
end

param_dir = [pwd '/params/' date_dir];
if ~exist(param_dir,'dir')
    mkdir(param_dir);
end

command_dir = [pwd '/commands/' date_dir];
if ~exist(command_dir,'dir')
    mkdir(command_dir);
end

rpinput_output_name = input_struct.sim_name;
rpinput_output_file = [rpinput_dir 'rpinput_' rpinput_output_name];

write = 1;
% check to see if you want to overwrite file
if exist(rpinput_output_file,'file')
   reply = input(['File ' rpinput_output_file ' exists. \n Do you want to overwrite? y/n '], 's');
   if (strcmp(reply,'n'))
      disp('Ok. That''s cool.');
      write = 0;
   end
end


% RPINPUT CALCULATOR
input_struct = CALC_RP(input_struct);

% RPINPUT WRITER
if write
    WRITE_RP(rpinput_template_file, rpinput_output_file, input_struct);
    save([param_dir 'param_' rpinput_output_name '.mat'], 'input_struct');
    run_dir = WRITE_CMD(command_dir, rpinput_output_name, input_struct.comp.mem,...
        input_struct.comp.tasks, input_struct.comp.run_time);
end

%exit;
