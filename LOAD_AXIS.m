function [x_axis, z_axis] = LOAD_AXIS(data_dir,data_type,file_number)

num_str = num2str(file_number,'%04d');
file    = [data_dir data_type '/' data_type '_' num_str '.h5'];

x_axis = h5read(file,'/AXIS/AXIS1'); % longitudinal axis start and stop in skin depths
z_axis = h5read(file,'/AXIS/AXIS2'); % transverse axis start and stop in skin depths