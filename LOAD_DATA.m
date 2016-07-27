function data = LOAD_DATA(data_dir,data_type,file_number)

if strncmp(data_type,'QEBI',4)
    num_str = '0000';
else
    num_str = num2str(file_number,'%04d');
end
    
loc     = ['/' data_type '-' num_str];
file    = [data_dir data_type '/' data_type '_' num_str '.h5'];

data = h5read(file,loc);