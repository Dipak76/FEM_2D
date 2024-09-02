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

% Clear the screen
clc;

% Clear all the variables from the workspace.
clear all ;

fclose('all');

%  The below file contains the instructions to set the input parameters
%  required for running the codes.

preprocessor

% The below file finds the elemental stiffness matrix, elemental external 
% force vector. Then it assembles the element matrices to get the global
% stifffness matrix and external force vector. Finally it solves for the
% displacement dU. Then it adds the displacement to undeformed coordinates
% Xn and obtains the deformed cooridnates xn.

solver

% file to save initial coordinates, deformed coocridnates, displacement and
% other specific outputs so desired in the Output/ folder.

fesave 

% The below file contains the function which plots the deformed and
% undeformed mesh. Also, on the defomred mesh you can also plot the
% contours of one of the following: 6 components of stress, 2 principal
% stresses, 3 invariants of stress, von Mises equivalent stress.
% The contour of displacement is plotted for dofs ux and uy .

 postprocessor






