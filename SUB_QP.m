clear all;

d = pwd;

MAKE_RP;

user = getenv('USER');
run_dir = ['/u/scratch/' user(1) '/' user '/' rpinput_output_name];

stat = mkdir(run_dir);
cd(run_dir);
copyfile(['/u/home/mori/sgess/executables/QuickPIC/' input_struct.comp.qpic_executable],'qpic.e');
copyfile([d '/rpinputs/' date_dir 'rpinput_' rpinput_output_name],'rpinput');
copyfile([d '/params/' date_dir 'param_' rpinput_output_name '.mat'],'param.mat');
copyfile([d '/commands/' date_dir 'qpic.e.cmd_' rpinput_output_name],'qpic.e.cmd');

system('qsub qpic.e.cmd','-echo');

cd(d);
exit;
