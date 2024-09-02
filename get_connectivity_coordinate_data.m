
% Element Connectivity
%
% Importing coordinates and connectivity matrix from mesh data exported from GMSH.

CON = [];

fid = fopen('Input/connectivity.txt','r') ;
while feof(fid) == 0
     CON = [CON;fscanf(fid,'%d\t %d\t %d\t %d',[1 4])];
end
fclose(fid);

% Initial Node Coordinates
Xn = [] ; % 2 dof per node i.e. x coordinate and y coordinate.

fid = fopen('Input/coordinate.txt','r') ;
while feof(fid) == 0
     Xn = [Xn;fscanf(fid,'%f\t %f',[1 2])];   
end
fclose(fid);

% Number of Elements
nel = length(CON) ;
% Number of Nodes
nno = length(Xn) ;

fid = fopen('Input/fe_data.txt','r') ;
FE = fscanf(fid,'%g');
ndoel = FE(1) ; % Number of dofs per element 
ndof = 2*nno ; ndoelo = ndoel ; % Number of dofs in the FE mesh, Number of dofs per element

fclose(fid);
clear FE;
% Current node coordinates initialized as initial coordinates to start
% with.
xn = Xn ; % Xn contains the coordinates of the reference configuration all the time
          % xn contains the coordinates of the current configuration all
          % the time