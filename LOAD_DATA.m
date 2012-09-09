function data = LOAD_DATA(data_dir,data_type,file_number)

num_str = num2str(file_number,'%04d');
loc     = ['/' data_type '-' num_str];
file    = [data_dir data_type '/' data_type '_' num_str '.h5'];

data = h5read(file,loc);