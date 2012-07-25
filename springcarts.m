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

==== ANALYSIS OF A SPRING-CART SYSTEM ====
This code performs three different types of analysis on the spring cart 
system developed during Lecture 1 of the FEA component of the course:
  1. A transient dynamic analysis.
  2. A static analysis.
  3. An undamped, unforced modal frequency analysis.
For the development of the system equations, please see the lecture slides.
The code demonstrates an application of the Direct Stiffness Method to a
particular example system.

== BOUNDARY CONDITIONS ==
= TRANSIENT DYNAMICS AND STATIC ANALYSIS =
For these analyses, the displacement of cart 1 is fixed at zero.  A force 
is applied in the positive direction to cart 4.
= MODAL FREQUENCY ANALYSIS =
In this analysis, no boundary conditions or forces are applied.

== UNITS ==
This example uses the following units:
  mm - length
  s  - time
  N  - force
  T  - mass

== OCTAVE NOTES ==
This code uses the ode45 function, which is not part of the main 
distribution of Octave.  To execute it, you will need the odepkg Octave
extension:
  http://octave.sourceforge.net/odepkg/
Instructions for installing Octave extensions on different platforms are 
available here:
  http://octave.sourceforge.net/

== MATLAB NOTES ==
Tested with Matlab R2012a, OS X.

%}


function springcarts()



  % If we're using Octave, set the graphics toolkit to the native 
  %  FLTK / OpenGL toolkit.  If this doesn't work, try commenting out the 
  %  graphics_toolkit() line.
  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  if isOctave
    graphics_toolkit('fltk');
  end


  %----
  % Compute the system matrices.
  %
  % This function simply assembles the system matrices in the form 
  % provided in the lecture.
  %
  % Inputs:
  %   k1, ..., k5 - spring stiffnesses
  %   c6          - viscous damping coefficient
  %   m1, ..., m4 - cart masses
  % Outputs:
  %   K - global system stiffness matrix
  %   C - global system damping matrix
  %   M - system mass matrix
  function [K,C,M] = systemMatrices(k1, k2, k3, k4, k5, c6, m1, m2, m3, m4)
    % Global system stiffness matrix
    K = [       k1          -k1           0            0
               -k1      k1+k2+k3+k4    -k2-k3         -k4
                0         -k2-k3      k2+k3+k5        -k5
                0           -k4          -k5         k4+k5     ];
    % Global system damping matrix
    Cz = zeros(4);
    Cz(1:2,1:2) = [ c6 -c6 ; -c6 c6 ];
    C = Cz;
    % System mass matrix
    M = diag([m1, m2, m3, m4]);
  end


  %----
  % Analysis of transient dynamics.
  %
  % This function performs forward integration of the system equations 
  % over a given time interval.  Fixed (time-invariant) boundary 
  % conditions are supported by enforcing zero velocity and acceleration 
  % on any non-moving degrees of freedom.
  %
  % Inputs:
  %   K        - global stiffness matrix (NxN) (N is the number of DOF)
  %   C        - global damping matrix (NxN)
  %   M        - global mass matrix (NxN)
  %   R        - vector of external forces (Nx1)
  %   tspan    - vector specifying the time interval of integration 
  %               [t0, tf]
  %   d0       - vector of initial displacements (Nx1)
  %   fixedDOF - boolean vector specifying any fixed degrees of freedom 
  %               (Nx1)
  % Outputs:
  %   T - vector of times at which the ODE has been solved
  %   d - displacement vector at each time sample in T
  function [T,d] = transientDynamics(K, C, M, R, tspan, d0, fixedDOF)

    % find the size of our system (nb: we are assuming that C and M are 
    %  the correct sizes, and are not doing any checks on them)
    [nDOF, n] = size(K);
    assert(n == nDOF, 'K must be a square matrix');

    % zero and identity matrices of the appropriate sizes for assembly
    Z = zeros(nDOF);
    I = eye(nDOF);

    % Mbar, Kbar and Rbar
    Mbar = [ I Z ; Z M ];
    Kbar = [ Z I ; -K -C ];
    Rbar = [ zeros(nDOF,1) ; R ];

    % Initial conditions (displacements are specified, velocities are 
    %  zero)
    q0 = [ zeros(nDOF,1) ; d0 ];

    % Matrix to set the derivatives of any fixed DOF to zero.  Mobility is 
    %  a matrix which multiplies the output derivatives of all "active" 
    %  DOF by 1, and all "fixed" or "inactive" DOF by 0.  If the first and
    %  second derivatives of a DOF are zero, then it will remain fixed 
    %  during the forward integration of the system equations, thus 
    %  enforcing the fixed boundary condition.
    MobilityD = 1 - (fixedDOF ~= 0);
    Mobility = diag([ MobilityD ; MobilityD ]);

    % First-order differential equation function for integration
    %  nb: this is an "anonymous function" because Octave does not support 
    %  nested functions
    qdot = @(t,q) Mobility * (Mbar \ (Kbar * q + Rbar));

    % Solve the ODE (nb: see the note in the comments at the start about 
    %  odepkg if you see an error here in Octave)
    options = odeset('RelTol', 1E-4, 'AbsTol', 1E-4, ...
                     'InitialStep', 0.1, 'MaxStep', 0.2);
    [T,q] = ode45(qdot, tspan, q0, options);

    % Pick out the values of q corresponding to the displacements d
    d = q(:,1:nDOF);

  end


  %-------------------------------%
  % Construct the system matrices %
  %-------------------------------%
  [K,C,M] = systemMatrices( ...
    ... % spring constants, k1...k6 (N/mm):
    5, 0.1, 0.1, 3, 0.2, ...
    ... % viscous damping coefficient (N/mm/s)
    1, ...
    ... % cart masses, m1...m4 (T)
    5E-3, 21E-3, 12E-3, 22E-3);


  %--------------------%
  % Transient dynamics %
  %--------------------%
  % Our external force is applied along DOF 4 (N):
  R = [ 0 0 0 100 ]';
  % The initial conditions for the DOFs are all zero (zero initial 
  %  displacement)
  d0 = [ 0 0 0 0 ]';
  % Our first DOF should be constant (fixed boundary condition)
  fixedDOF = [ 1 0 0 0 ]';
  % Let's set the final time through which we should integrate the 
  %  equations (s)
  tFinal = 10;
  % Perform the integration
  [T,dTrans] = transientDynamics(K, C, M, R, [0 tFinal], d0, fixedDOF);


  %-----------------%
  % Static analysis %
  %-----------------%
  % Create d_E, the "essential" boundary condition displacements, for 
  %  which the nodal displacements are known.
  %  In this case, d_E = [ u1 ] = [ 0 ]
  d_E = 0;
  % Now create the force vector for the free nodes.  In this 
  %  case, f_F = [ R2 R3 R4 ] = [ 0 0 R4 ]
  f_F = R(2:end);
  % Next, partition K into [ KE KEF ; KEF' KF ].  We're taking a special 
  %  case here, since we know that we have one essential boundary
  %  condition corresponding to the first DOF.  This would need to be more 
  %  general in the case of other constraints.
  KE  = K(1,1);
  KEF = K(1,2:end);     % row 1, columns 2->end
  KF  = K(2:end,2:end); % rows 2->end, columns 2->end
  % Solve for the unknown static displacements
  d_F = KF \ (f_F - (KEF') * d_E);
  % Solve for the unknown reaction at node 1
  rE_1 = KE * d_E + KEF * d_F;


  %--------------------------------%
  % Plot dynamic vs static results %
  %--------------------------------%
  figure();
  hold on;
  % plot the dynamics response curves
  plot(T,dTrans(:,1), 'b');
  plot(T,dTrans(:,2), 'g');
  plot(T,dTrans(:,3), 'r');
  plot(T,dTrans(:,4), 'm');
  % plot the static response lines
  plot([0 tFinal], [d_F(1), d_F(1)], 'g--');
  plot([0 tFinal], [d_F(2), d_F(2)], 'r--');
  plot([0 tFinal], [d_F(3), d_F(3)], 'm--');
  legend('u1', 'u2', 'u3', 'u4');
  xlabel('Time (s)');
  ylabel('Displacement (mm)');
  hold off;


  %-------------------%
  % Modal frequencies %
  %-------------------%
  % Solve the eigenvalue problem
  [d_Modal, lambda] = eig(M\K);
  % Compute the angular frequencies (the square root of each lambda value)
  omega = diag(lambda) .^ (0.5);
  % Convert to cycles per time frequency
  freq = sort(real(omega ./ (2*pi)));
  % Display results
  disp('Un-damped, un-forced modal frequencies (Hz):');
  disp(freq);


  %-----------------------------------%
  % Export data (for carts animation) %
  %-----------------------------------%
  %fps = 30
  %subframes = 16
  %tResample = [0:(1/(30*16)):tFinal];
  %dResample = interp1(T, dTrans, tResample);
  %dlmwrite('carts.csv', dResample, ',');


end
