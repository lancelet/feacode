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

function truss2d()

  function [Ke] = trussStiffness2D(x1, y1, x2, y2, A, E)
    dx = x2 - x1;
    dy = y2 - y1;
    le = sqrt(dx^2 + dy^2);   % element length
    ke = A * E / le;          % element spring constant
    phie = atan2(dy, dx);     % element angle
    % small rotation matrix
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

  function [g] = computeGlobalNode(x, n1, n2)
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
      eData = elements(e,1:end);
      n1 = eData(1); n2 = eData(2); A = eData(3); E = eData(4);
      Ke = trussStiffness2D( ...
        nodes(n1,1), nodes(n1,2), nodes(n2,1), nodes(n2,2), A, E);
      for i=1:4,
        gi = computeGlobalNode(i, n1, n2);
        li = i;
        for j=1:4,
          gj = computeGlobalNode(j, n1, n2);
          lj = j;
          K(gi, gj) = K(gi, gj) + Ke(li, lj);
        end
      end
    end
  end 

  function plotTruss(nodes, elements, d, F, dMag, fMag)
    % number of elements
    elementSize = size(elements);
    nElements = elementSize(1);
    % number of ndoes
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
      gx1 = computeGlobalNode(1, n1, n2);
      gy1 = computeGlobalNode(2, n1, n2);
      gx2 = computeGlobalNode(3, n1, n2);
      gy2 = computeGlobalNode(4, n1, n2);
      % compute strains, stresses and forces from original displacements
      qi = nodes(n2,1:end) - nodes(n1,1:end);
      qf = qi + [ d(gx2) d(gy2) ] - [ d(gx1) d(gy1) ];
      le2 = qi * qi';
      lf2 = qf * qf';
      stress = (sqrt(lf2 / le2) - 1) * E;
      force = stress * A;
      fprintf(1, 'Element %2d; F = % 6.02f kN, stress = % 6.02f MPa\n', ...
        e, force / 1E3, stress / 1E6);
      % magnified displacements (for plotting)
      x1 = nodes(n1,1) + dMag * d(gx1);
      y1 = nodes(n1,2) + dMag * d(gy1);
      x2 = nodes(n2,1) + dMag * d(gx2);
      y2 = nodes(n2,2) + dMag * d(gy2);
      plot([x1, x2], [y1, y2], 'k', 'LineWidth', 1.5, 'LineSmoothing', 'on');
    end
    for n=1:nNodes,
      gx = computeGlobalNode(1, n, n);
      gy = computeGlobalNode(2, n, n);
      x = nodes(n, 1) + dMag * d(gx);
      y = nodes(n, 2) + dMag * d(gy);
      fx = fMag * F(gx);
      fy = fMag * F(gy);
      quiver([x x], [y y], [fx 0], [0 fy], ...
        'Color', [0.8, 0, 0], 'LineWidth', 1.0, 'AutoScale', 'off', ...
        'MaxHeadSize', 5);
    end
    xlim([-3 10]);
    ylim([-4, 6]);
    set(gca, 'DataAspectRatio', [1 1 1]);
    hold off;
  end

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

  fixedDOF = 3;

  K = directAssembly(nodes, elements);
  K
  d_E = [ zeros(fixedDOF,1) ];
  f_F = [ zeros(13 - 3,1) ; -14000 ];
  KE = K(1:fixedDOF,1:fixedDOF);
  KEF = K(1:fixedDOF,(fixedDOF+1):end);
  KF = K((fixedDOF+1):end,(fixedDOF+1):end);
  d_F = KF \ (f_F - (KEF') * d_E);
  r_E = KE * d_E + KEF * d_F;
  d = [ d_E ; d_F ];
  F = [ r_E ; f_F ];

  plotTruss(nodes, elements, d, F, 200, 0.00005);

end
