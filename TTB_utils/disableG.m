function disableG(h,k)

    if get(h(k),'Value')==2
        set(h(k+3),'String',num2str(0.5));
        set(h(k+3),'Visible','off');
        set(h(k+4),'String',num2str(0.5));
        set(h(k+4),'Visible','off');
        
        %msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
        %get(h(k))
        
    else
    
        get(h(k),'Value')~=2
        set(h(k+3),'String',num2str(0));
        set(h(k+3),'Visible','on');
        set(h(k+4),'String',num2str(3));
        set(h(k+4),'Visible','on');
        
        %msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
        %get(h(k))
        
    end

end