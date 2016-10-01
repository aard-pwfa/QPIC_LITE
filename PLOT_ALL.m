function list = PLOT_ALL(data_dir,data_dim,file_number,param_struct,plot_dir,save_plot,save_ext,plot_some)

cmap = custom_cmap;
x = dir(data_dir);
list = cell(0,1);

for i = 1:numel(x)
    if ~x(i).isdir; continue; end;
    if strncmp(x(i).name,'.',1); continue; end;
    if strcmp(x(i).name,'LOG'); continue; end;
    if strcmp(x(i).name,'RAW-BEAM'); continue; end;
    if x(i).name(4)=='I'; continue; end;
    list{end+1} = x(i).name;
    
end

if nargin == 8
    use = false(size(list));
    for i=1:numel(list)
        if sum(strcmp(list{i},plot_some))
            use(i) = true;
        end
    end
    list = list(use);
end

for i = 1:numel(list)
    
    type = list{i};
    
    data = LOAD_DATA(data_dir, type, file_number);
    PLOT_FUN(type,data,data_dim,data_dir,cmap.bwr,file_number,param_struct,plot_dir,save_plot,save_ext,i);
    
end

if sum(strcmp('QEB-XZ',list))
    typeB = 'QEB-XZ';
    typeP = 'QEP1-XZ';
    typeT = 'QET-XZ';
    dataB = LOAD_DATA(data_dir, typeB, file_number);
    dataP = -LOAD_DATA(data_dir, typeP, file_number);
    [~,b] = max(abs(dataB(:)));
    if dataB(b) > 0; dataB = -dataB; end;
    data = 10*dataB+dataP;
    PLOT_FUN(typeT,data,data_dim,data_dir,cmap.bwr,file_number,param_struct,plot_dir,save_plot,save_ext,i+1);
end
if sum(strcmp('QEB-YZ',list))
    typeB = 'QEB-YZ';
    typeP = 'QEP1-YZ';
    typeT = 'QET-YZ';
    dataB = LOAD_DATA(data_dir, typeB, file_number);
    dataP = -LOAD_DATA(data_dir, typeP, file_number);
    [~,b] = max(abs(dataB(:)));
    if dataB(b) > 0; dataB = -dataB; end;
    data = dataB+dataP;
    PLOT_FUN(typeT,data,data_dim,data_dir,cmap.bwr,file_number,param_struct,plot_dir,save_plot,save_ext,i+2);
end
if sum(strcmp('QEB-XY',list))
    typeB = 'QEB-XY';
    typeP = 'QEP1-XY';
    typeT = 'QET-XY';
    dataB = LOAD_DATA(data_dir, typeB, file_number);
    dataP = -LOAD_DATA(data_dir, typeP, file_number);
    [~,b] = max(abs(dataB(:)));
    if dataB(b) > 0; dataB = -dataB; end;
    data = dataB+dataP;
    PLOT_FUN(typeT,data,data_dim,data_dir,cmap.bwr,file_number,param_struct,plot_dir,save_plot,save_ext,i+3);
end

if sum(strcmp('FEX-XZ',list))
    typeE = 'FEX-XZ';
    typeB = 'FBY-XZ';
    typeW = 'FWX-XZ';
    dataE = LOAD_DATA(data_dir, typeE, file_number);
    dataB = LOAD_DATA(data_dir, typeB, file_number);
    
    data = dataE-dataB;
    PLOT_FUN(typeW,data,data_dim,data_dir,cmap.bwr,file_number,param_struct,plot_dir,save_plot,save_ext,i+4);
end
if sum(strcmp('FEX-YZ',list))
    typeE = 'FEX-YZ';
    typeB = 'FBY-YZ';
    typeW = 'FWX-YZ';
    dataE = LOAD_DATA(data_dir, typeE, file_number);
    dataB = LOAD_DATA(data_dir, typeB, file_number);
    
    data = dataE-dataB;
    PLOT_FUN(typeW,data,data_dim,data_dir,cmap.bwr,file_number,param_struct,plot_dir,save_plot,save_ext,i+5);
end
if sum(strcmp('FEX-XY',list))
    typeE = 'FEX-XY';
    typeB = 'FBY-XY';
    typeW = 'FWX-XY';
    dataE = LOAD_DATA(data_dir, typeE, file_number);
    dataB = LOAD_DATA(data_dir, typeB, file_number);
    
    data = dataE-dataB;
    PLOT_FUN(typeW,data,data_dim,data_dir,cmap.bwr,file_number,param_struct,plot_dir,save_plot,save_ext,i+6);
end