function [rel] = compute_distributed_load(ndoel,con,xigs,ngps,xI,boundary,t)

rel = zeros(ndoel,1) ;

ie = [con(1)-1 con(1) con(2)-1 con(2) con(3)-1 con(3) con(4)-1 con(4)] ;

if isempty(intersect(ie,boundary))==0 % if intersection of current element's nodes and boundary nodes is not zero
    for gp=1:ngps
        xi=xigs(gp,1);
        wg=xigs(gp,2);

        N1=(1-xi)/2;
        N2=(1+xi)/2;
        N=[N1 N2];

        dpN=[-1/2 1/2];
        J1=0;
        J2=0;
        for i=1:2
            J1 = J1 + dpN(1,i)*xI(4-i,1); % for top boundary, jacobian to transform x to xi
            J2 = J2 + dpN(1,i)*xI(1+i,2); % for right boundary, jacobian to transform y to eta
        end
        
        if isempty( intersect( ie,boundary(:,2) ) )==0
            rel(3,1) = rel(3,1) + N(1)*t(2)*J2*wg;   % assigning force along x dof to element's southeast node
            rel(5,1) = rel(5,1) + N(2)*t(2)*J2*wg;   % assigning force along x dof to element's northeast node
        end
        if isempty( intersect( ie,boundary(:,1) ) )==0
            rel(6,1) = rel(6,1) + N(1)*t(1)*J1*wg;  % assigning force along y dof to element's northwest node
            rel(4,1) = rel(4,1) + N(2)*t(1)*J1*wg;  % assigning force along y dof to element's northeast node
        end
    end
 end  

