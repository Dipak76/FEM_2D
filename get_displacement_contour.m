for j=[1,2]
pu = j; % 1 = displacement along x ; 2 = displacement along y
figure (j+4);

    for i = 1:nel
                %{
                 Xp = [ PX(CON(i,1),1) PX(CON(i,2),1) ; PX(CON(i,4),1) PX(CON(i,3),1) ] ;
                 Yp = [ PX(CON(i,1),2) PX(CON(i,2),2) ; PX(CON(i,4),2) PX(CON(i,3),2) ] ; 
                 Cp = [ Ux(CON(i,1),pu)  Ux(CON(i,2),pu) ;  
                          Ux(CON(i,4),pu)  Ux(CON(i,3),pu) ] ;  

                 s=pcolor(Xp,Yp,Cp) ; hold on
                 set(s,'EdgeColor','none');
                %}
              Xn1 = PX(CON(i,1),:) ;
              Xn2 = PX(CON(i,2),:) ;   
              Xn3 = PX(CON(i,3),:) ;   
              Xn4 = PX(CON(i,4),:) ;
              Xn5 = ( Xn1 + Xn2 + Xn3 + Xn4 )/4 ;

              Cp1 = Ux(CON(i,1),pu) ;
              Cp2 = Ux(CON(i,2),pu) ;
              Cp3 = Ux(CON(i,3),pu) ;
              Cp4 = Ux(CON(i,4),pu) ;
              Cp5 = ( Cp1 + Cp2 + Cp3 + Cp4)/4 ;

              Xp = [ Xn1(1)  Xn2(1) ; Xn5(1)  Xn5(1) ] ;
              Yp = [ Xn1(2)  Xn2(2) ; Xn5(2)  Xn5(2) ] ;
              Cp = [ Cp1     Cp2    ; Cp5     Cp5    ] ;
              s=pcolor(Xp,Yp,Cp) ; hold on ;
              set(s,'EdgeColor','none');

              Xp = [ Xn2(1)  Xn3(1) ; Xn5(1)  Xn5(1) ] ;
              Yp = [ Xn2(2)  Xn3(2) ; Xn5(2)  Xn5(2) ] ;
              Cp = [ Cp2     Cp3    ; Cp5     Cp5    ] ;
              s=pcolor(Xp,Yp,Cp) ; hold on ;
              set(s,'EdgeColor','none');

              Xp = [ Xn3(1)  Xn4(1) ; Xn5(1)  Xn5(1) ] ;
              Yp = [ Xn3(2)  Xn4(2) ; Xn5(2)  Xn5(2) ] ;
              Cp = [ Cp3     Cp4    ; Cp5     Cp5    ] ;
              s=pcolor(Xp,Yp,Cp) ; hold on ;
              set(s,'EdgeColor','none');

              Xp = [ Xn4(1)  Xn1(1) ; Xn5(1)  Xn5(1) ] ;
              Yp = [ Xn4(2)  Xn1(2) ; Xn5(2)  Xn5(2) ] ;
              Cp = [ Cp4     Cp1    ; Cp5     Cp5    ] ;
              s=pcolor(Xp,Yp,Cp) ; hold on ;
              set(s,'EdgeColor','none');
    end

    % Figure Properties
    xlabel('X / L_o')
    ylabel('Y / L_o')
    if j==1
        string_title = strcat('x-Deformation Contour');
    elseif   j==2
        string_title = strcat('y-Deformation Contour');
    end
    title(string_title) ;

    colorbar
    axis(axe)
    axis equal
end