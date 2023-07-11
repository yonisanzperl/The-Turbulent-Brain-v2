    function disablePerts(h,k)

    if get(h(k),'Value')
        
            
        set(h(k+1),'String',{'errHete','Turbu','Transfer','IC'});
        set(h(k+1),'Visible','on');
        set(h(k+2),'String',num2str(3));
        set(h(k+2),'Visible','on');
        
        
        
        %msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
        %get(h(k))
        
    else

        set(h(k+1),'String',num2str(0.5));
        set(h(k+1),'Visible','off');
        set(h(k+2),'String',num2str(0.5));
        set(h(k+2),'Visible','off');
        set(h(k+3),'String',num2str(0.5));
        set(h(k+3),'Visible','off');
         set(h(k+4),'String',num2str(0.5));
        set(h(k+4),'Visible','off');
         set(h(k+5),'String',num2str(0.5));
        set(h(k+5),'Visible','off');
        %msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
        %get(h(k))
        
    end

end