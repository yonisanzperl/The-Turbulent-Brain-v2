    function disablePerts(h,k)

    if get(h(k),'Value')
        
            
        set(h(k+2),'String',num2str(0.18));
        set(h(k+2),'Visible','on');
        set(h(k+3),'String',num2str(10));
        set(h(k+3),'Visible','on');
        
        
        %msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
        %get(h(k))
        
    else

        set(h(k+2),'String',num2str(0.5));
        set(h(k+2),'Visible','off');
        set(h(k+3),'String',num2str(0.5));
        set(h(k+3),'Visible','off');
 
        
    end

end