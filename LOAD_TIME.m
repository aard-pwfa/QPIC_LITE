function [iter, dt, time] = LOAD_TIME(data_dir,data_type,file_number)

num_str = num2str(file_number,'%04d');
file    = [data_dir data_type '/' data_type '_' num_str '.h5'];

iter = h5readatt(file,'/','ITER');
dt   = h5readatt(file,'/','DT');
time = h5readatt(file,'/','TIME');