% PLot the stress contour
% Current Stresses

festrs11

figure(10);  %plot original mesh(undeformed)
for i = 1:nel
    Xdev = [Xn(CON(i,1),1) Xn(CON(i,2),1) Xn(CON(i,3),1) Xn(CON(i,4),1) Xn(CON(i,1),1)] ;
    Ydev = [Xn(CON(i,1),2) Xn(CON(i,2),2) Xn(CON(i,3),2) Xn(CON(i,4),2) Xn(CON(i,1),2)] ;        
    
    plot(Xdev,Ydev,'k') ; hold on
    fill(Xdev,Ydev,'b') ;
end
xlabel('X / L_o')
ylabel('Y / L_o')
title('Undeformed Plate') ;

axis(axe)
axis equal

get_displacement_contour

% Compute the exact stresses of the infinite plate with 
% a circular centered hole problem
get_exact_stresses