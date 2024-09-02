%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 INDIAN INSTITUTE OF TECHNOLOGY GUWAHATI                 %
%                  DEPARTMENT OF MECHANICAL ENGINEERING                   %
%                                                                         %
%                               2023-24 3rd SEMESTER                           %
%                                                                         %
%                           FINITE ELEMENT METHODS                 %
%                                                                         %
%                                                                         %
% Code initially developed by: Sachin Singh Gautam                        %
%                                                                         %
%                                                                         %
% Assignment 2: Due date 31.07.2024, Friday, 5 PM                            %
%                                                                         %
% The code is written for finding the displacement, strains and % 
%  stresses for a rectangular plate with hole subjected to biaxial
%  tension applied in the form of uniformly distributed load    %
%  For finite element simulation the code used 4-noded quadrilateral %
%  element.             
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FELP Stress Recovery 


% Initialize
stress_at_nodes  = zeros(nno,6) ;
KL = zeros(nno,1) ;

% Element Loop
for i = 1:nel   

    % Previous Timestep Dofs
    con = 2*CON(i,:) ;
    ie = [con(1)-1 con(1) con(2)-1 con(2) con(3)-1 con(3) con(4)-1 con(4)] ;
    
    % Current Nodal Coordinates and local Displacements
    oI = Xn(CON(i,:),:) ;
    xI = xn(CON(i,:),:) ;
    uI = U(ie) ; 
    
    % Loop over Gauss Points 
    sel = zeros(4,4) ;
    bel = zeros(4,2) ;
    mll = zeros(4,1) ;
    
    
    % Gauss Point Loop (Volume Integrals)
    for gp = 1:ngpv
 
        % Gauss point coordinates, weight
        xi = xigv(gp,1) ; eta = xigv(gp,2) ; wg = xigv(gp,3) ;
        
        % Shape Functions
        N1 = ( 1 - xi )*( 1 - eta )/4 ;
        N2 = ( 1 + xi )*( 1 - eta )/4 ;
        N3 = ( 1 + xi )*( 1 + eta )/4 ;
        N4 = ( 1 - xi )*( 1 + eta )/4 ;
        N  = [N1 N2 N3 N4] ;
         
      % Shape Fct (xi,eta)-Derivatives
        dpN = [ -( 1-eta )  -( 1-xi ) ;
                 ( 1-eta )  -( 1+xi ) ;
                 ( 1+eta )   ( 1+xi ) ;
                -( 1+eta )   ( 1-xi ) ]/4 ;

        % Jacobian between Current Configuration and Master Element
        Jac = zeros(2,2) ;
        for j = 1:4
           Jac = Jac + xI(j,:)' * dpN(j,:) ;
        end
        detJ = det(Jac) ; 

        % Deformation Warning
        if detJ < 0       
            jterm = 1 ;
        end

        invJ = inv(Jac) ;

        % Shape Fct (x,y)- Spatial Derivatives
        for j = 1:4
           dN(j,:)  =  dpN(j,:) * invJ ;
        end
        
    % B-Matrix           
    B = [ dN(1,1)        0  dN(2,1)        0  dN(3,1)        0  dN(4,1)        0 ;
                0  dN(1,2)        0  dN(2,2)        0  dN(3,2)        0  dN(4,2) ;
          dN(1,2)  dN(1,1)  dN(2,2)  dN(2,1)  dN(3,2)  dN(3,1)  dN(4,2)  dN(4,1) ] ;
      
         % Deformation Gradient
    invF = zeros(2,2) ;
    for j = 1:4
       invF = invF + oI(j,:)' * dN(j,:) ; 
    end
    Fgr  = inv(invF) ; 
    detF = det(Fgr) ;
    
    if detF < 0 
        jterm = 1 ;              
    end  
    
        strain_at_gp = B * uI ;
        stress_at_gp = D * strain_at_gp ;
        sigma_33 = nu * (stress_at_gp(1,1) + stress_at_gp(2,1) ) ;      
        sig = [ stress_at_gp(1,1) ; stress_at_gp(2,1) ; sigma_33 ; stress_at_gp(3,1) ] ; 
        
        % Stress Projection
        sel = sel + N'*sig' * detJ * wg ;
        
        % Body Forces projected onto the deformed Configuration
        bel = bel + N' * ldv*g/detF * detJ * wg ;
        
        % 'Mass' Matrix
        mll = mll + N'   * detJ * wg ;
        
    end
    
    % Structure Arrays
    stress_at_nodes(CON(i,:),:)        = stress_at_nodes(CON(i,:),:)        + [sel bel] ; 
    KL(CON(i,:),1)       = KL(CON(i,:),1)       + mll ;
    
end

% Lumped Mass Solution
for i = 1:nno
    stress_at_nodes(i,:) = stress_at_nodes(i,:)/KL(i,1) ;
end

% Initialize
stress_at_nodes2  = zeros(nno,12) ;
stress_at_nodes2(:,1:6) = stress_at_nodes(:,:) ;

% Principal Stress
stress_at_nodes2(:,7)  = (stress_at_nodes2(:,1)+stress_at_nodes2(:,2))/2 + sqrt( ((stress_at_nodes2(:,1)-stress_at_nodes2(:,2)).^2)/4 + stress_at_nodes2(:,4).^2 ) ;
stress_at_nodes2(:,8)  = (stress_at_nodes2(:,1)+stress_at_nodes2(:,2))/2 - sqrt( ((stress_at_nodes2(:,1)-stress_at_nodes2(:,2)).^2)/4 + stress_at_nodes2(:,4).^2 ) ;
    
% Invariants
stress_at_nodes2(:,9)  = stress_at_nodes2(:,1) + stress_at_nodes2(:,2) + stress_at_nodes2(:,3) ;
stress_at_nodes2(:,10) = stress_at_nodes2(:,7).*stress_at_nodes2(:,8) + stress_at_nodes2(:,8).*stress_at_nodes2(:,3) + stress_at_nodes2(:,3).*stress_at_nodes2(:,7) ;
stress_at_nodes2(:,11) = stress_at_nodes2(:,7).*stress_at_nodes2(:,8).*stress_at_nodes2(:,3) ;
    
% Von Mises
stress_at_nodes2(:,12) = sqrt(( (stress_at_nodes2(:,7)-stress_at_nodes2(:,8)).^2 + (stress_at_nodes2(:,8)-stress_at_nodes2(:,3)).^2 + (stress_at_nodes2(:,3)-stress_at_nodes2(:,7)).^2 )/2) ;

% File for printing the deformed  coordinates
fileStressAtNodes = 'Output/stress_at_nodes.txt' ;
% Open Mesh File
fid = fopen(fileStressAtNodes,'w') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Current Configuration and Stresses
for i = 1:nno
    fprintf(fid,'%17.10e %17.10e %17.10e %17.10e %17.10e %17.10e  %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e  \n',...
    stress_at_nodes(i,1:6),stress_at_nodes2(i,7:12));
end

fclose(fid) ;

pstrs = 1 ;   % 1 - means sigma_xx; 2 - means sigma_yy; 3 - means sigma_zz
              % 4 - means sigma_xy; 5 - means sigma_yz; 6 - means sigma_zx
              % 7 and 8 - principal stresses
              % 9 - means 1st invariant of stress; 10 - means 2nd invariant of
              % stress; 11: means 3rd invariant of stress
              % 12 - means von Mises stress
              
% Load the data file
PX1 = load('Output/initial_coordinate.txt') ;
PX2 = load('Output/deformed_coordinate.txt') ;
PX3 = load('Output/deformed_displacement.txt') ;

stress_file_name = 'Output/stress_at_nodes.txt' ;
stress_at_nodes_plot = load(stress_file_name) ;
 
% FELP Deformed Mesh Plot
% stress_at_nodesoger Sauer - 11.8.08

jsmooth = 1 ;
PXDeformed = PX2 + factor * PX3 ;

% Coordinates to plot
PX = PXDeformed ;

% Plot Elements  
for j=[1,2,3,4,12]    
figure(j)
pstrs = j;
    for i = 1:nel

      if pstrs == 0   % Plot Colored Mesh

        X = [PX(CON(i,1),1) PX(CON(i,2),1) PX(CON(i,3),1) PX(CON(i,4),1) PX(CON(i,1),1)] ;
        Y = [PX(CON(i,1),2) PX(CON(i,2),2) PX(CON(i,3),2) PX(CON(i,4),2) PX(CON(i,1),2)] ;        

        if col == 'o'
            plot(X,Y,'k') ; hold on
        else
    %         if  con_array(i,1) == 1 
             fill(X,Y,col) ; hold on    
    %         else
    %            fill(X,Y,'r') ; hold on            
    %         end
        end

      else  % Plot Stresses onto Mesh

        % Four-noded, incorrect Interpolation  
        if jsmooth == 0 

          Xp = [ PX(CON(i,1),1) PX(CON(i,2),1) ; PX(CON(i,4),1) PX(CON(i,3),1) ] ;
          Yp = [ PX(CON(i,1),2) PX(CON(i,2),2) ; PX(CON(i,4),2) PX(CON(i,3),2) ] ;        
          Cp = [ stress_at_nodes_plot(CON(i,1),pstrs)  stress_at_nodes_plot(CON(i,2),pstrs) ;  
              stress_at_nodes_plot(CON(i,4),pstrs)  stress_at_nodes_plot(CON(i,3),pstrs) ] ;  

          s=pcolor(Xp,Yp,Cp) ; hold on ; %shading interp
          set(s,'EdgeColor','none');

        % Three-noded, smoother Interpolation  
        elseif jsmooth == 1

          Xn1 = PX(CON(i,1),:) ;
          Xn2 = PX(CON(i,2),:) ;   
          Xn3 = PX(CON(i,3),:) ;   
          Xn4 = PX(CON(i,4),:) ;
          Xn5 = ( Xn1 + Xn2 + Xn3 + Xn4 )/4 ;

          Cp1 = stress_at_nodes_plot(CON(i,1),pstrs) ;
          Cp2 = stress_at_nodes_plot(CON(i,2),pstrs) ;
          Cp3 = stress_at_nodes_plot(CON(i,3),pstrs) ;
          Cp4 = stress_at_nodes_plot(CON(i,4),pstrs) ;
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

        % Highly smooth Intrpolation (expensive)  
        else  
          Conp = CON' ; plsig = [0 pstrs] ; 
          fesmip2

        end
      end
    end

    MaxLx = max(PX2(:,1));
    MaxLy = max(PX2(:,2));
    Mx = 0.1*Lx ; axe = [Xo-Mx Xo+MaxLx+Mx Yo-Mx Yo+MaxLy+Mx] ;

    % Figure Properties
    xlabel('X / L_o')
    ylabel('Y / L_o')
    colorbar
    if j==1
        string_title = strcat('\sigma_{xx} distribution (Numerical solution)');
        caxis( [0 3.5e7] );
    elseif j==2
        string_title = strcat('\sigma_{yy} distribution (Numerical solution) ');
        caxis( [0 10e7] );
    elseif j==3
        string_title = strcat('\sigma_{zz} distribution ');
        caxis( [0 3.1e7] );
    elseif j==4
        string_title = strcat('\sigma_{xy} distribution {Numerical solution} ');
        caxis( [-5e7 5e7] );
    elseif j==12
        string_title = strcat('Equivalent(Von-Mises) Stress distribution - Numerical solution ');
        caxis( [0 4.5e7] );
    end
    title(string_title) ;

    %shading interp
    axis(axe)
    axis equal
end