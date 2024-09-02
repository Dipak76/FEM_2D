
% --------> 1
% File for printing the initial coordinates
fileInitialCoordinate = 'Output/initial_coordinate.txt' ;
 % Opening the file for the first time
fidInitialCoordinate = fopen(fileInitialCoordinate,'w') ;
for i = 1:nno % run the loop over all the nodes
    fprintf(fidInitialCoordinate,'%20.15f \t %20.15f \n',Xn(i,1),Xn(i,2));
end
fclose(fidInitialCoordinate) ; % Close the file

% --------> 2
% File for printing the deformed  coordinates
fileDeformedCoordinate = 'Output/deformed_coordinate.txt' ;
% Opening the file for the first time
fidDeformedCoordinate = fopen(fileDeformedCoordinate,'w') ; 
for i = 1:nno
    fprintf(fidDeformedCoordinate,'%20.15f \t %20.15f \n',xn(i,1),xn(i,2));
end
fclose(fidDeformedCoordinate) ; % Close the file


% --------> 3
% File for printing the deformed  displacement
fileDeformedDisplacement = 'Output/deformed_displacement.txt' ;
% Opening the file for the first time
fidDeformedDisplacement = fopen(fileDeformedDisplacement,'w') ; 
for i = 1:nno % run the loop over all the nodes
    fprintf(fidDeformedDisplacement,'%20.15f \t %20.15f \n',Ux(i,1),Ux(i,2));
end
fclose(fidDeformedDisplacement) ; % Close the file
