function beam = LOAD_BEAM(data_dir,file_number)

num_str = num2str(file_number,'%04d');
file    = [data_dir '/RAW-BEAM/01/RAW-BEAM-01_' num_str '.h5'];

x = h5read(file,'/x1');
y = h5read(file,'/x2');
z = h5read(file,'/x3');

px = h5read(file,'/p1');
py = h5read(file,'/p2');
pz = h5read(file,'/p3');

beam = [x y z px py pz];