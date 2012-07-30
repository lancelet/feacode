%-------------------------------------------------------------------------
% MCEN40024 SOLID MECHANICS (2012)
%  DEPARTMENT OF MECHANICAL ENGINEERING
%  THE UNIVERSITY OF MELBOURNE
%     ___________ 
%    / __/ __/ _ |
%   / _// _// __ |
%  /_/ /___/_/ |_|  - Finite Element Analysis
%
% This code may be used freely under the Creative Commons Attribution 
% Unported (CC BY 3.0) license:  
%   http://creativecommons.org/licenses/by/3.0/
% Cite the attribution as: 
%  "Portions of this code written by Dr Jonathan Merritt"
%
% Original Author: Jonathan Merritt <j.s.merritt@gmail.com>
%                                   <merritt@unimelb.edu.au>
%-------------------------------------------------------------------------

%{
  
==== ANALYSIS OF 2D TRUSSES ====

(aka: First year Mech Eng student in 507 lines of code...)

This code performs static analyses of 2D trusses, with elements that 
should be equivalent to T2D2 elements in Abaqus.  It is set up to solve 
the truss problem provided in Lecture 2 of the FEA component of the 
course.  It also solves problem #1.3.32 from the Abaqus 6.8
verification manual, which tests 2D truss elements.

The truss is described below by the `nodes' and `elements' variables, 
which specify the geometry of the problem.  Then a global stiffness 
matrix corresponding to the specified geometry is assembled using the
direct assembly technique and the element equations for 2D trusses that
were described during the lecture.  The static system equation is 
solved using the partitioning method.

==== BOUNDARY CONDITION LIMITATION ====

The code can only apply fixed (zero displacement) boundary conditions
to individual DOF.  Any DOF that are fixed MUST APPEAR FIRST in the 
list of global DOF.  The global DOF are arranged as follows:

  Global DOF # | DOF
             1   node 1 x1
             2   node 1 x2
             3   node 2 x1
             4   node 2 y1
            ...   ...

The fixedDOF variable below specifies the number of DOF that are fixed.
In this case, fixedDOF = 3, which means:

  Global DOF # | DOF       | Fixed or Free?
             1   node 1 x1 | Fixed
             2   node 1 x2 | Fixed
             3   node 2 x1 | Fixed
             4   node 2 y1 | Free
             5   node 3 x1 | Free
            ...   ...         ...

This limitation arises because of the matrix partitioning that is used
for the solution.  It could be circumvented by separate assembly of
the matrices required for the solution.

== MATLAB NOTES ==
Tested with Matlab R2012a, OS X.

== OCTAVE NOTES ==
Untested.

%}

function truss2d()


  %----
  % Compute the local / compacted stiffness matrix for a 2D truss element.
  %
  % This follows the derivation provided in the lecture.
  %
  % Inputs:
  %   x1, y1, x2, y2 - coordinates of the ends of the 2D truss element
  %   A              - area of the truss element
  %   E              - Young's modulus of the element
  % Outputs:
  %   Ke - local / compacted stiffness matrix for the element
  function [Ke] = trussStiffness2D(x1, y1, x2, y2, A, E)

    % properties of the element
    dx = x2 - x1;             % length in x
    dy = y2 - y1;             % length in y
    le = sqrt(dx^2 + dy^2);   % element length
    ke = A * E / le;          % element spring constant
    phie = atan2(dy, dx);     % element angle

    % small rotation matrix (converts back to an x-aligned truss element)
    Rbare = [
       cos(phie) sin(phie)
      -sin(phie) cos(phie)
    ];
    % block rotation matrix
    Re = [ Rbare zeros(2,2) ; zeros(2,2) Rbare ];
    
    % stiffness matrix for un-rotated element
    Kprime = [ 
       1 0 -1 0
       0 0  0 0
      -1 0  1 0
       0 0  0 0
    ];
    
    % stiffness matrix for rotated element
    Ke = ke * Re' * Kprime * Re;

  end


  %----
  % Computes local to global mapping of DOF for our truss problem.
  %
  % NB: THIS FUNCTION ONLY WORKS FOR ELEMENTS WITH 4 DOF EACH!
  %
  % This function is a mapping between local DOF of an element and the global
  % DOF.  Suppose there is an element which spans global nodes 3 and 5; for
  % this element:
  %   n1 = 3
  %   n2 = 5
  % The element itself has 4 DOF (x=1 through x=4).  In the 2D truss elements, 
  % the DOF correspond to x1, y1, x2, y2.  The global DOF corresponding to 
  % these local DOF are:
  %   x = 1 (x1) :  global node # = 2*n1-1 = 2*3-1 = 5
  %   x = 2 (y1) :  global node # = 2*n1   = 2*3   = 6
  %   x = 3 (x2) :  global node # = 2*n2-1 = 2*5-1 = 9
  %   x = 4 (y2) :  global node # = 2*n2   = 2*5   = 10
  % So it's a simple packing where two DOF are stored for each global node
  % number used to reference them.  The table for all DOF looks like this:
  %  DOF # | DOF
  %      1   node 1 x1
  %      2   node 1 y1
  %      3   node 2 x1
  %      4   node 2 y1
  %     ...   ...
  %   2N-1   node N x1
  %     2N   node N y1
  %
  % Inputs:
  %   x  - the _local_ DOF of the element (1, 2, 3, 4)
  %   n1 - the number of the first global node of the element
  %   n2 - the number of the second global node of the element
  % Outputs:
  %   g  - the _global_ DOF corresponding to _local_ DOF x
  function [g] = computeGlobalDOF(x, n1, n2)
    switch x
    case 1
      g = 2 * n1 - 1;
    case 2
      g = 2 * n1;
    case 3
      g = 2 * n2 - 1;
    case 4
      g = 2 * n2;
    end
  end


  %-----
  % Performs direct assembly of the system stiffness matrix.
  %
  % Inputs:
  %   nodes - A matrix containing nodes of the form:
  %       [ x1 y1 ; x2 y2 ; x3 y3 ; ... ]
  %     where (x1,y1) are the coordinates of the first node, (x2,y2) are the
  %     coordinates of the second, etc.
  %   elements - A matrix containing 2D truss elements of the form:
  %       [ n1^1 n2^2 A^1 E^1 ; n1^2 n2^2 A^2 E^2 ; ... ]
  %     where n1 and n2 are 1-based indexes into the array of nodes, and
  %     specify the global connectivity of the element.  A is the cross-
  %     sectional area and E is the Young's modulus
  %
  % Outputs:
  %   K - the system stiffness matrix, computed by direct assembly
  function [K] = directAssembly(nodes, elements)

    % number of elements
    elementSize = size(elements);
    nElements = elementSize(1);
    % 2D problem; there are 2 x nodes DOF
    nodeSize = size(nodes);
    nDOF = nodeSize(1) * 2;

    % set initial K
    K = zeros(nDOF, nDOF);

    % iterate over elements
    for e=1:nElements,
    
      % fetch data from the array of elements
      eData = elements(e,1:end);
      n1 = eData(1); n2 = eData(2); A = eData(3); E = eData(4);
    
      % compute the local / compacted stiffness matrix for the element
      Ke = trussStiffness2D( ...
        nodes(n1,1), nodes(n1,2), nodes(n2,1), nodes(n2,2), A, E);
    
      % scatter the compacted stiffness matrix using direct assembly
      for i=1:4,    % iterate over rows of K
        gi = computeGlobalDOF(i, n1, n2);   % global DOF
        li = i;                             % local DOF
        for j=1:4,  % iterate over columns of K
          gj = computeGlobalDOF(j, n1, n2); % global DOF
          lj = j;                           % local DOF
          K(gi, gj) = K(gi, gj) + Ke(li, lj);
        end
      end
    
    end
  end 


  %----
  % Plots a truss.
  %
  % Plots undeformed and deformed shapes, and external forces.
  %
  % Inputs:
  %   nodes    - matrix containing nodes (see the directAssembly function)
  %   elements - matrix containing elements (see the directAssembly function)
  %   d        - column matrix of global node displacements
  %   F        - column matrix of global node external forces
  %   dMag     - magnification factor for plotting the displacements
  %   fMag     - magnification factor for plotting the forces
  function plotTruss(nodes, elements, d, F, dMag, fMag)

    % number of elements
    elementSize = size(elements);
    nElements = elementSize(1);
    % number of nodes
    nodeSize = size(nodes);
    nNodes = nodeSize(1);
    
    % plot the original (un-deformed) truss
    figure();
    hold on;
    for e=1:nElements,
      eData = elements(e,1:end);
      n1 = eData(1); n2 = eData(2); A = eData(3); E = eData(4);
      plot([nodes(n1,1), nodes(n2,1)], [nodes(n1,2), nodes(n2,2)], ...
        'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.5, 'LineSmoothing', 'on');
    end

    % plot the deformed truss and forces
    for e=1:nElements,
      eData = elements(e,1:end);
      n1 = eData(1); n2 = eData(2); A = eData(3); E = eData(4);
      gx1 = computeGlobalDOF(1, n1, n2);
      gy1 = computeGlobalDOF(2, n1, n2);
      gx2 = computeGlobalDOF(3, n1, n2);
      gy2 = computeGlobalDOF(4, n1, n2);
      % magnified displacements (for plotting)
      x1 = nodes(n1,1) + dMag * d(gx1);
      y1 = nodes(n1,2) + dMag * d(gy1);
      x2 = nodes(n2,1) + dMag * d(gx2);
      y2 = nodes(n2,2) + dMag * d(gy2);
      plot([x1, x2], [y1, y2], 'k', 'LineWidth', 1.5, 'LineSmoothing', 'on');
    end
    
    % plot the forces acting on the nodes
    for n=1:nNodes,
      gx = computeGlobalDOF(1, n, n);
      gy = computeGlobalDOF(2, n, n);
      x = nodes(n, 1) + dMag * d(gx);
      y = nodes(n, 2) + dMag * d(gy);
      fx = fMag * F(gx);
      fy = fMag * F(gy);
      quiver([x x], [y y], [fx 0], [0 fy], ...
        'Color', [0.8, 0, 0], 'LineWidth', 1.0, 'AutoScale', 'off', ...
        'MaxHeadSize', 5);
    end

    % fix up the plot dimensions
    xlim([-3 10]);
    ylim([-4, 6]);
    set(gca, 'DataAspectRatio', [1 1 1]);
    hold off;

  end


  %----
  % Inspects a truss; printing out forces and stresses.
  %
  % To compute the internal stresses and forces in the truss, it is tempting
  % simply to compute the strain as (change in length) / (length) and use the 
  % material constitutive laws to find the stress and force.  However, that's 
  % the wrong way to go about things.
  %
  % To compute the internal stresses and forces, the same constitutive 
  % relations and interpolations that applied to the original element should
  % be used.  To achieve that here, the original element stiffness matrices
  % are re-computed, and the relationship F = Kd is used with the computed
  % displacements to find F.  F and A are then used to find stress.
  %
  % Inputs:
  %   nodes    - matrix containing nodes (see the directAssembly function)
  %   elements - matrix containing elements (see the directAssembly function)
  %   d        - column matrix of global node displacements
  %   F        - column matrix of global node external forces
  function inspectTruss(nodes, elements, d, F)

    % number of elements
    elementSize = size(elements);
    nElements = elementSize(1);
    % number of ndoes
    nodeSize = size(nodes);
    nNodes = nodeSize(1);

    % print out forces and stresses
    fprintf('Force and stress (-ve is compressive)\n');
    for e=1:nElements,
    
      % extract the element data
      eData = elements(e,1:end);
      n1 = eData(1); n2 = eData(2); A = eData(3); E = eData(4);
    
      % find which global DOF correspond to our local DOF x1, y1, x2, y2
      gx1 = computeGlobalDOF(1, n1, n2);
      gy1 = computeGlobalDOF(2, n1, n2);
      gx2 = computeGlobalDOF(3, n1, n2);
      gy2 = computeGlobalDOF(4, n1, n2);
    
      % re-compute the stiffness matrix for the element
      x1 = nodes(n1,1); y1 = nodes(n1,2); x2 = nodes(n2,1); y2 = nodes(n2,2);
      Ke = trussStiffness2D(x1, y1, x2, y2, A, E);
    
      % use F = Kd to find the nodal forces
      Fe = Ke * [ d(gx1) ; d(gy1) ; d(gx2) ; d(gy2) ];
    
      % dot-product the nodal force against a normal vector along the original
      %  element direction (bit hacky, but it works)
      F1 = [ Fe(1) Fe(2) ];
      dl = [ (x2-x1) (y2-y1) ]';
      nhat = dl / sqrt(dl' * dl);
      Fmag = -F1 * nhat;
    
      % put force in kN and stress in MPa, and print them
      F_kn = Fmag / 1E3; % kN
      stress_kPa = F_kn / A / 1E3; % MPa
      fprintf(1, 'e:%2d | F = % 8.4f kN | stress = % 10.6f MPa \n', ...
        e, F_kn, stress_kPa);
    
    end
  end

  %----
  % Computes a static FEA solution.
  %
  % The system equation F=Kd is partitioned as follows:
  %    [ r_E   = [  KE  KEF   [ d_E
  %      f_F ]     KEF'  KF ]   d_F ]
  % Where 
  %   d_E - known (specified) displacement BCs (Nx1)
  %   d_F - unknown displacement BCs (Mx1)
  %   r_E - unknown nodal forces (Nx1)
  %   f_F - known (specified) nodal forces (Mx1)
  %    KE - (NxN), part of stiffness matrix K
  %   KEF - (NxM), part of stiffness matrix K
  %    KF - (MxM), part of stiffness matrix K
  %
  % The second row of the matrix is solved first to find the nodal 
  % displacements:
  %   d_F = inv(KF) * (f_F - KEF' * d_E)
  % Then the first row is used to find the unknown external forces
  %   r_E = KE * d_E + KEF * d_F
  %
  % Inputs:
  %   d_E - known nodal displacements (essential BCs) (Nx1)
  %   f_F - known nodal forces (natural BCs) (Mx1)
  %   K   - system stiffness matrix ((N+M)x(N+M))
  % Outputs:
  %   d - nodal displacements ((N+M)x1)
  %   F - nodal forces ((N+M)x1)
  function [d,F] = solveStatic(d_E, f_F, K)

    % get the size of the important matrices and DOF
    f_F_size = size(f_F);
    d_E_size = size(d_E);
    K_size = size(K);
    nDOF = K_size(1);
    nEssentialDOF = d_E_size(1);

    % check that some important sizes are all compatible
    assert(K_size(2) == nDOF, 'K must be square');
    assert(f_F_size(1) == (nDOF - nEssentialDOF), ...
      '# of free DOF must be (nDOF - nEssentialDOF): problem with f_F');

    % compute the solution
    KE = K(1:nEssentialDOF,1:nEssentialDOF);
    KEF = K(1:nEssentialDOF,(nEssentialDOF+1):end);
    KF = K((nEssentialDOF+1):end,(nEssentialDOF+1):end);
    d_F = KF \ (f_F - (KEF') * d_E);
    r_E = KE * d_E + KEF * d_F;

    % Re-assemble, Stephanie
    d = [ d_E ; d_F ];  % global displacement
    F = [ r_E ; f_F ];  % global external node forces

  end


  %---- SOLUTION 1 : TRUSS FROM LECTURE 2 ----

  fprintf(1, '---- Solution 1 : Truss from lecture 2 ----\n');

  % nodes contains the (x,y) values of nodal coordinates, one per row in (m)
  nodes = [ ...
    0 3
    0 0
    3 3
    3 0
    6 3
    6 1
    9 3
  ];

  % some element properties
  A = pi * 0.02 ^ 2;  % area of all truss elements, in m^2
  E = 200E9;          % Young's modulus of all truss elements, in Pa

  % elements contains the parameters of each element:
  %  node 1, node 2, area, Young's modulus
  elements = [ ...
    1 2 A E
    2 4 A E
    1 4 A E
    1 3 A E
    3 4 A E
    4 6 A E
    3 6 A E
    3 5 A E
    5 6 A E
    5 7 A E
    6 7 A E
  ];

  % how many total DOF are there?
  nodeSize = size(nodes);
  nDOF = 2 * nodeSize(1);

  % how many fixed DOF in the problem?  for this code, all the fixed DOF have
  %  to appear first in the list of DOF (see comments at the top)
  fixedDOF = 3;

  % use direct assembly to find the system stiffness matrix
  K = directAssembly(nodes, elements);

  % solve the system
  d_E = [ zeros(fixedDOF,1) ];  % essential boundary conditions
  f_F = [ zeros(nDOF - fixedDOF - 1,1) ; -14000 ];  % applied forces
  [d,F] = solveStatic(d_E, f_F, K);

  % make pretty pictures and dump solution information
  plotTruss(nodes, elements, d, F, 200, 0.00005);
  inspectTruss(nodes, elements, d, F);


  %---- SOLUTION 2 : ABAQUS VERIFICATION PROBLEM 1.3.32 FOR T2D2 ELEMENTS ----
  % Variables have a _v appended, for `verification'

  fprintf(1, '\n\n');
  fprintf(1, '---- Solution 2 : Abaqus verification #1.3.32 ----\n');

  fprintf(1, 'Reference Solution (from the Abaqus verification manual):\n');
  fprintf(1, 'e: 1 | stress = 0.032907 MPa\n');
  fprintf(1, 'e: 2 | stress = 0.041134 MPa\n');
  fprintf(1, 'e: 3 | stress = 0.032907 MPa\n');

  nodes_v = [ ...
     5 10
    10 10
    15 10
    10  0
  ];  

  A_v = 0.1;
  E_v = 30E6;
  elements_v = [ ...
    1 4 A_v E_v
    2 4 A_v E_v
    3 4 A_v E_v
  ];

  nodeSize_v = size(nodes_v);
  nDOF_v = 2 * nodeSize_v(1);  
  fixedDOF_v = 6;

  K_v = directAssembly(nodes_v, elements_v);

  d_E_v = [ zeros(fixedDOF_v,1) ]; % essential boundary conditions
  f_F_v = [ zeros(nDOF_v - fixedDOF_v - 1,1) ; -10000 ];  % applied forces
  [d_v,F_v] = solveStatic(d_E_v, f_F_v, K_v);

  inspectTruss(nodes_v, elements_v, d_v, F_v);


end
