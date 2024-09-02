% Normalization Parameters

fid = fopen('Input/normalization_parameters.txt','r') ;
NP = fscanf(fid,'%f \t %f \t %f');

Lo = NP(1); % Length
Eo = NP(2); % Youngs Modulous
To = NP(3); % Time

fclose(fid);
clear NP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find the origin of the problem

fid = fopen('Input/origin.txt','r') ;
OR = fscanf(fid,'%f \t %f');

Xo = OR(1) ; 
Yo = OR(2) ;

fclose(fid);
clear OR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Dimension of the geometry; 
% 
fid = fopen('Input/geometric_data.txt','r') ;
DIM = fscanf(fid,'%f \t %f  \t %f \t %f');

Lx = DIM(1); % Lx = Length of the plate in the x direction
Ly = DIM(2); % Ly = Length of the plate in y direction
R = DIM(3); % R = Radius of hole
thickness_of_plate = DIM(4); % thickness_of_plate = this is thickness of the plate in the z direction.Note that this has to be changed for plane stress cases 

fclose(fid);
clear DIM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Material properties of the plate
%
fid = fopen('Input/material_data.txt','r') ;
MP = fscanf(fid,'%f \t %f \t %g');

E = MP(1) ; % E = Young's modulus of the plate
nu = MP(2) ; % nu = Poissions ratio of the plate
option_type_2D = MP(3) ;% option_type_2D = 1 --> plane strain % option_type_2D = 2 --> plane stress

fclose(fid);
clear MP;

get_material_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finite element data

get_connectivity_coordinate_data

% Displacement rearragement arrays - this is needed at the end when the
% displacement is added to the current coordinates to get the deformed
% coordinates

i1 = zeros(nno,1) ; i2 = zeros(nno,1) ;
for i = 1:nno
    i1(i) = 2*i-1 ; % x dof of all the nodes
    i2(i) = 2*i ; % y dofs of all the nodes
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get the boundary condition array
%
get_boundary_conditon_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% External Force Vector - this has to be computed for distributed load
F_ext = zeros(ndof,1) ;  

% for distributed load
top_boundary=[];
right_boundary=[];
for i=1:nno
    if Xn(i,2) == Ly
        top_boundary=[top_boundary;2*i];
    end
    if Xn(i,1) == Lx
        right_boundary=[right_boundary;2*i-1];
    end
end
boundary=[top_boundary right_boundary];

fid = fopen('Input/external_loads.txt','r') ;
EL = fscanf(fid,'%f \t %f \t %g \t %f \t %f ');

traction_top = EL(1) ; % E = Young's modulus of the plate
traction_right = EL(2) ; % nu = Poissions ratio of the plate
ldv = EL(3) ; % Flag to enable gravity load
g = [EL(4) EL(5)]; % acceleration values

fclose(fid);
clear EL;

t=[traction_top traction_right];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial displacements. U is total displacement whereas DU is incremental
% displacement
U = zeros(ndof,1) ; DU = zeros(ndof,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadature parameters for volume intergation
jgaus = 2 ; fequad ; xigv = xig2 ; ngpv = jgaus^2 ; 

% Quadature parameters for surface intergation. This is needed when the
% surface integrals in cases like UDL are required.
jgaus = 5 ; fequad ; xigs = xig1 ; ngps = jgaus ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Call the function to find the size of arrays needed for parallel
% computing so that code can used for meshes with extremely large number of
% elements.

parallel_computing_array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Factor by which the displacement must be multiplied to get a visible
% deformed geometry at the end of the deformation process

fid = fopen('Input/factor.txt','r') ;
factor = fscanf(fid,'%g');
fclose(fid);

 % FLag for error catching. 0 means no error. It will be set to nonzero 
 % value in case an error is detected
jterm = 0 ;

% Computation Options
jinfo = 2 ;
jcomp = 1 ;