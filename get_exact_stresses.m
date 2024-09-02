exact_stresses = zeros(nno,3);

for i =1:nno  % computing exact stresses
    
    r = sqrt(xn(i,1)*xn(i,1) + xn(i,2)*xn(i,2));
    theta = atan(xn(i,2)/xn(i,1));
    c2t = cos(2*theta);
    c4t = cos(4*theta);
    s2t = sin(2*theta);
    s4t = sin(4*theta);

    fac1 = (R/r)^2;
    fac2 = fac1*fac1;

    exact_stresses(i,1) = t(2) - 0.5*fac1*c2t*( 3*t(2) - t(1) ) + ( t(2) - t(1) )*( 1.5*fac2 - fac1 )*c4t;  % here t(2) is traction on right boundary t(1) is traction on top boundary
    exact_stresses(i,2) = t(1) - 0.5*fac1*c2t*( t(2) - 3*t(1) ) - ( t(2) - t(1) )*( 1.5*fac2 - fac1 )*c4t;
    exact_stresses(i,3) = ( t(1) + t(2) )*( (1.5*fac2-fac1)*s4t ) - ( t(2) - t(1) )*(0.5*fac1*s2t);  
    
end    
for j=1:3
    pstre=j;

    figure(j+6);   % plot exact stresses on deformed mesh
    for i=1:nel
        %{
        Xp = [ PX(CON(i,1),1) PX(CON(i,2),1) ; PX(CON(i,4),1) PX(CON(i,3),1) ] ;
        Yp = [ PX(CON(i,1),2) PX(CON(i,2),2) ; PX(CON(i,4),2) PX(CON(i,3),2) ] ;        
        Cp = [  exact_stresses(CON(i,1),pstre)   exact_stresses(CON(i,2),pstre) ;  
                   exact_stresses(CON(i,4),pstre)   exact_stresses(CON(i,3),pstre) ] ;  

        s=pcolor(Xp,Yp,Cp) ; hold on ; %shading interp
        set(s,'EdgeColor','none');
        %}
          Xn1 = PX(CON(i,1),:) ;
          Xn2 = PX(CON(i,2),:) ;   
          Xn3 = PX(CON(i,3),:) ;   
          Xn4 = PX(CON(i,4),:) ;
          Xn5 = ( Xn1 + Xn2 + Xn3 + Xn4 )/4 ;

          Cp1 = exact_stresses(CON(i,1),pstre) ;
          Cp2 = exact_stresses(CON(i,2),pstre) ;
          Cp3 = exact_stresses(CON(i,3),pstre) ;
          Cp4 = exact_stresses(CON(i,4),pstre) ;
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
    colorbar
    if j==1
        string_title = strcat('\sigma_{xx} distribution (Analytical Solution) ');
        caxis( [0 3.5e7] );
    elseif j==2
        string_title = strcat('\sigma_{yy} distribution (Analytical Solution) ');
        caxis( [0 10e7] );
    else
        string_title = strcat('\sigma_{xy} distribution (Analytical Solution) ');
        caxis( [-5e7 5e7] );
    end
    title(string_title) ;

   %shading interp
    axis(axe)
    axis equal
end