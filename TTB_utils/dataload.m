function dataload(h,k)

%load(get(handles(k),'String'))

%set(handles(k+2),'items', {'a';'b'})

% 


[datos_name,path_emp] = uigetfile('*.mat','Elegir los datos empíricos'); %carga datos empíricos

pathfinal_datos=[path_emp datos_name];

datos_listvar = who( matfile(pathfinal_datos));

%if size(datos_listvar,1)==1
    %datos_listvar{end+1,1}='No hay nada cargado';
%end


set(h(k+1),'Value',1);

set(h(k+1),'String',datos_listvar );

set(h(k+2),'String',pathfinal_datos);




end