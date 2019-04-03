function [M,K,DOF,Omega,Phi] = FE_fun(x)
%
% Parametrized FE-assembly of a 2D Offshore wind turbine on a jacket support structure
% including a simplified foundation model 
%
% Required completions: 
% Assign the elements of "x" to a parameter of the FE-model in the "Assign input" section:
% 	To overwrite the desired parameters you can change "***" to the names of the desired
%	variables from the code below and comment the corresponding variables in the code itself. 
%
% Function input: 
% - x: vector containing variables for optimisation 
%  
% Function Output: 
% - M: element mass matrix
% - K: element stiffness matrix
% - DOF: dof list
% - Omega: list of natural frequencies [rad/s]
% - Phi: matrix containing the eigenmodes
%
% Notes: all visual output surpressed
%
% Authors: Pim van der Male & Dominik Fallais
% -------------------------------------------------------------------------

%% Assign input 
% assign input variables to selected parameters
%Esoil = x(1);     	% change *** to a variable from code below 
%Econ = x(2);     	% and comment that variable (by adding a "%" in front of the line) 

var1 = evalin('base','var1');
var2 = evalin('base','var2');

%% Geometry definition
%    -                 ####       --> turbine mass (rotor nacelle assembly)
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |htow              ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    |                  ||
%    -                ======         --> transition piece
%    |                |\  /|         
%    |                | \/ |         --> segment 4
%    |                | /\ |             (leg 2) 
%    |                |/  \|           
%    |               |\    /|          
%    |               |  \/  |        --> segment 3
%    |               |  /\  |            (leg 2) 
%    |htot           |/    \|          
%    |              |\      /|         
%    |              |  \  /  |       --> segment 2         
%    |              |  /  \  |           (leg 2) 
%    |              |/      \|         
%    |  -          |\        /|        
%    |  |          |   \  /   |      --> segment 1        
%    |  |hseg      |   /  \   |          (leg 1) 
%    |  |          |/        \|        
%    -  -      x---+----------+---x  --> bottom member and horizontal soil springs     
%				   |		  |      --> vertical soil spring					
%				   x		  x		 --> boundary conditions: x = clamped.						
%                      btop
%                     |----|
%                  |----------|
%                      bbot 
%
%   Notes on Foundation moddelling:
%   The soil-structure interaction is modelled by using 1D spring (rod)  
%   elements in the vertical and horizontal direction at each pile tip.
%   Since no accurate model could be set up for the foundation (due to lack
%   of information), a simple model is implemented; the parameters of this 
%   model do not represent real soil properties which could be measured in 
%   the field, but equivalent/effective properties. Instead of having a 
%   distributed stiffness and mass, the rigidity and mass of the soil and 
%   pile are captured in equivalent "lumped" parameters which are assigned 
%   to a single spring in each direction; for simplicity the vertical and 
%   horizontal springs are assumed to be equal. 
%
%   Soil spring parameters:  (geometry): Afnd; (material): Esoil, nusoil, 
%   rhosoil. 

% Foundation
% ----------
Afnd = [2.5; 2.5];    % effective Area factor for simplified spring foundation [vertical, horizontal] 
nfnd = 2;             % number of foundation sections

% Turbine
% -------

% Turbine mass definition
Mnac = 240000;                                                                                                   % nacelle mass [kg]
Mra = 110000;                                                                                                    % rotor assembly mass [kg]
Mtop = Mnac + Mra;                                                                                               % turbine mass [kg]

switch var1
    case 4
        Mtop = x(1);
end
switch var2
    case 4
        Mtop = x(2);
end

% Jacket 
% ------

% Global geometry of the jacket support structure 
htot = 70;                                                                                                       % jacket height [m]
nseg = 4;                                                                                                        % brace segment number [-]
hseg = htot / nseg;                                                                                              % brace segment height [m]
bbot = 12;                                                                                                       % jacket base width [m]
btop = 8;                                                                                                        % jacket top width [m]

% Jacket width distribution over the height
b = zeros(nseg + 1,1);

for j = 1:nseg + 1
    b(j,1) = bbot - (j - 1) * (bbot - btop) * hseg / htot;                                                       % jacket width at top segment j [m]
end

% Tubular sections geometry
% Three different tubular sections are distinguished for the definition of the jacket support structure.
nts = 3;                                                                                                         % number of different tubular sections [-]
D(1,1) = 1.289;                                                                                                  % diameter leg 1 [m]
D(2,1) = 1.123;                                                                                                  % diameter leg 2 [m]
D(3,1) = 0.732;                                                                                                  % diameter brace / bottom member [m]
tw(1,1) = 0.0537;                                                                                                % wall thickness leg 1 [m]
tw(2,1) = 0.0312;                                                                                                % wall thickness leg 2 [m]
tw(3,1) = 0.020;                                                                                                 % wall thickness brace / bottom member [m]

Ats = zeros(nts,1);                         
Its = zeros(nts,1);

% The 2D definition of the jacket support structure requires the definition of combined leg and brace sections. 
Ats(1,1) = 2 * 0.25 * pi * (D(1,1)^2 - (D(1,1) - 2 * tw(1,1))^2 + D(3,1)^2 - (D(3,1) - 2 * tw(3,1))^2);          % combined deep leg + brace cross-sectional area [m^2]
Ats(2,1) = 2 * 0.25 * pi * (D(2,1)^2 - (D(2,1) - 2 * tw(2,1))^2 + D(3,1)^2 - (D(3,1) - 2 * tw(3,1))^2);          % combined leg + brace cross-sectional area [m^2]
Ats(3,1) = 2 * 0.25 * pi * (D(3,1)^2 - (D(3,1) - 2 * tw(3,1))^2);                                                % combined brace cross-sectional area [m^2]
Its(1,1) = 2 * 0.25 * pi * ((0.5 * D(1,1))^4 - (0.5 * (D(1,1) - 2 * tw(1,1)))^4);                                % combined deep leg moment of inertia [m^4]
Its(2,1) = 2 * 0.25 * pi * ((0.5 * D(2,1))^4 - (0.5 * (D(2,1) - 2 * tw(2,1)))^4);                                % combined leg moment of inertia [m^4]
Its(3,1) = 2 * 0.25 * pi * ((0.5 * D(3,1))^4 - (0.5 * (D(3,1) - 2 * tw(3,1)))^4);                                % combined brace moment of inertia [m^4]

% Transition piece
% ----------------

% Transition piece dimensions
btp = 9.6;                                                                                                       % transition piece width [m]
htp = 5.25;                                                                                                      % transition piece height [m]
Atp = btp * htp;                                                                                                 % transition piece cross-sectional area [m^2]
Itp = 1/12 * btp * htp^3;                                                                                        % transition piece moment of inertia [m^4]

% Tower
% -----

% Global geometry of the turbine tower
htow = 70.5;                                                                                                     % tower height [m]
nstow = 9;                                                                                                       % tower segment number [-]
hstow = htow / nstow;                                                                                            % tower segment height [m]
Dbot = 5.5;                                                                                                      % tower bottom diameter [m]
tbot = 0.034;                                                                                                    % tower bottom wall thickness [m]
Dtop = 4.0;                                                                                                      % tower top diameter [m]
ttop = 0.020;                                                                                                    % tower top wall thickness [m]

% Tower diameter and wall thickness distribution over the height
Dtow = zeros(nstow,1);
ttow = zeros(nstow,1);

for j = 1:nstow
    D(j+nts,1) = Dbot - (j - 1) * (Dbot - Dtop) * hstow / htow;                                                  % tower diameter at segment j [m]
    tw(j+nts,1) = tbot - (j - 1) * (tbot - ttop) * hstow / htow;                                                 % tower wall thickness at segment j [m]
end

Atow = zeros(nstow,1);
Itow = zeros(nstow,1);

for j = 1:nstow
    Atow(j,1) = 0.25 * pi * (D(j+nts,1)^2 - (D(j+nts,1) - 2 * tw(j+nts,1))^2);                                   % tower cross-sectional area at segment j [m^2]
    Itow(j,1) = 0.25 * pi * ((0.5 * D(j+nts,1))^4 - (0.5 * (D(j+nts,1) - 2 * (tw(j+nts,1))))^4);                 % tower moment of inertia at segment j [m^4]
end    



%% Section Defintion 
% Sections matrix
% ---------------

% Sections = [SecID A    ky   kz   Ixx Iyy Izz]
% - SecID: section ID [-]
% - A: section cross-sectional area [m^2]
% - ky: shear deflection factor [-]
% - kz: shear deflection factor [-]
% - Ixx: moment of inertia [m^4]
% - Iyy: moment of inertia [m^4]
% - Izz: moment of inertia [m^4]

Sections = zeros(nfnd + nts + nstow + 1,7);

% Tubular sections
for j = 1:nts
    Sections(j,1) = j;
    Sections(j,2) = Ats(j,1);
    Sections(j,3) = Inf;
    Sections(j,4) = Inf;
    Sections(j,5) = 0;
    Sections(j,6) = 0;
    Sections(j,7) = Its(j,1);
end

% Transition piece
Sections(nts + 1,1) = nts + 1;
Sections(nts + 1,2) = Atp;
Sections(nts + 1,3) = Inf;
Sections(nts + 1,4) = Inf;
Sections(nts + 1,5) = 0;
Sections(nts + 1,6) = 0;
Sections(nts + 1,7) = Itp;

% Tower
for j = 1:nstow
    Sections(nts +1 + j,1) = nts + 1 + j;
    Sections(nts +1 + j,2) = Atow(j,1);
    Sections(nts +1 + j,3) = Inf;
    Sections(nts +1 + j,4) = Inf;
    Sections(nts +1 + j,5) = 0;
    Sections(nts +1 + j,6) = 0;
    Sections(nts +1 + j,7) = Itow(j,1);
end


% Foundation 
% ----------
sfnd = [max(find(Sections(:,1)))+1  max(find(Sections(:,1)))+2];      % segment numbers for foundation == last segment number +1 // +2

for j = 1:nfnd
    Sections(nts +1 + nstow + j,:) = [sfnd(j)  Afnd(j,1) NaN NaN NaN NaN NaN];
end

%% Material properties

% steel properties
% ----------------
Est = 210e9;                                                                                                        % Young's modulus steel [N/m^2]
nust = 0.3;                                                                                                         % poisson coeff. steel 
rhost = 7850;                                                                                                       % density steel [kg/m^3]

% concrete properties
% -------------------
Econ = 10e9;                                                                                                        % Young's modulus concrete [N/m^2]
nucon = 0.2;                                                                                                        % poisson coeff. concrete
rhocon = 2485;                                                                                                      % density concrete [kg/m^3]
                                                                                                                    
% soil properties
% ---------------
Esoil = 3e9;                                                                                                    % A Young's modulus of soil [N/m^2] 
nusoil = 0.1;                                                                                                       % A poisson coeff. soil
rhosoil = 1250;                                                                                                     % A Density of soil [kg/m^3]

switch var1
    case 1
        Esoil = x(1);
    case 2
        Est = x(1);
    case 3
        rhost = x(1);
end
switch var2
    case 1
        Esoil = x(2);
    case 2
        Est = x(2);
    case 3
        rhost = x(2);
end

% Materials matrix
% ----------------
% Materials = [MatID E nu rho]
% - MatID: material ID [-]
% - E: Young's modulus [N/m^2]
% - nu: Poisson's ratio [-]
% - rho: mass density [kg/m^3]

Materials = [1 Est   nust   rhost;
             2 Econ  nucon  rhocon; 
             3 Esoil nusoil rhosoil];  

%% FE segment division
nNfnd = 2;                                                                                                         % FE division soil springs: 2 pairs of springs [-] 
nNbot = 1;                                                                                                         % FE division bottom member [-]
nNtop = 4;                                                                                                         % FE division top member [-]
nNseg = 2;                                                                                                         % FE division leg segment [-]
nNtow = 9;                                                                                                         % FE division tower [-]

% Element types matrix
% --------------------
% Types = {EltTypID    EltName}
% - EltTypID: element type ID [-]
% - EltName: element type [-] 

Types = {1 'beam';
         2 'truss'};

%% Node definition

% Node matrix
% -----------

% Nodes = [NID X Y Z]
% - NID: node ID
% - X: node X coordinate [m]
% - Y: node Y coordinate [m]
% - Z: node Z coordinate [m]

Nodes = zeros(500,4);                                                                                               % preallocate nodes
count = 1;                                                                                                          % counter representing node ID
 

% Botom member node definition
% ----------------------------

Nbm = zeros(nNbot+1,1);                                                                                             % botom member node ID array
for k = 0:nNbot
    Nodes(count,1) = count;
    Nodes(count,2) = -b(1,1) / 2 + k * b(1,1) / nNbot;
    Nodes(count,3) = 0;
    Nodes(count,4) = 0;
    Nbm(k+1,1) = count;
    count = count + 1;
end

% Legs and braces node definition
% -------------------------------

Nleg1 = zeros(nNseg,nseg);                                                                                          % leg node ID array
Nleg2 = zeros(nNseg,nseg);                                                                                          % leg node ID array
Nbr1 = zeros(nNseg + 1,nseg);                                                                                       % brace node ID array
Nbr2 = zeros(nNseg + 1,nseg);                                                                                       % brace node ID array
cmid = zeros(nseg,1);                                                                                               % brace midpoint node ID array

for j = 1:nseg  % Leg node definition
    
    if j == 1
       Nleg1(1,j) = Nodes(1,1);
       Nbr1(1,j) = Nleg1(1,j);
       Nleg2(1,j) = Nodes(nNbot + 1,1);
       Nbr2(1,j) = Nleg2(1,j);
    else
       Nleg1(1,j) = Nleg1(1,j-1)+nseg;
       Nbr1(1,j) = Nbr1(1,j-1)+nseg;
       Nleg2(1,j) = Nleg2(1,j-1)+nseg;
       Nbr2(1,j) = Nbr2(1,j-1)+nseg;
    end
    
    for k = 1:nNseg
        Nodes(count,1) = count;
        Nodes(count,2) = -b(j,1) / 2 + k * (b(j,1)-b(j + 1,1)) / 2 / nNseg;
        Nodes(count,3) = 0;
        Nodes(count,4) = (j - 1) * hseg + k * hseg / nNseg;
        Nleg1(k + 1 ,j) = count;
        Nbr2(3,j) = Nleg1(k + 1 ,j);
        count = count + 1;

        Nodes(count,1) = count;
        Nodes(count,2) = b(j,1) / 2 - k * (b(j,1)-b(j + 1,1)) / 2 / nNseg;
        Nodes(count,3) = 0;
        Nodes(count,4) = (j - 1) * hseg + k * hseg / nNseg;
        Nleg2(k + 1 ,j) = count;
        Nbr1(3,j) = Nleg2(k + 1 ,j);
        count = count + 1;
    end
end
for j = 1:nseg  % Brace node definition
    
    cmid(j,1) = count;
    Nodes(count,1) = cmid(j,1);
    Nodes(count,2) = 0;
    Nodes(count,3) = 0;
    Nodes(count,4) = (j - 1) * hseg + b(j,1) / (b(j + 1,1) + (b(j,1) - b(j + 1,1)) / 2) * hseg / 2;
    count = count + 1;

    Nbr1(2,j) = cmid(j,1);
    Nbr2(2,j) = cmid(j,1);
end
     

% Transition piece node definition
% --------------------------------

Ntp = zeros(nNtop + 3,1);                                                                                           % transition piece node ID array

Nodes(count,1) = count;
Nodes(count,2) = -btp / 2;
Nodes(count,3) = 0;
Nodes(count,4) = htot;
Ntp(1,1) = count;
count = count + 1;

Ntp(2,1) = max(Nleg1(:,4));
Ntp(nNtop + 2,1) = max(Nleg2(:,4));

Nodes(count,1) = count;
Nodes(count,2) = btp / 2;
Nodes(count,3) = 0;
Nodes(count,4) = htot;
Ntp(nNtop + 3,1) = count;
count = count + 1;

for k = 1:nNtop - 1
    Nodes(count,1) = count;
    Nodes(count,2) = -b(5,1) / 2 + k * b(5,1) / nNtop;
    Nodes(count,3) = 0;
    Nodes(count,4) = htot;
    Ntp(k + 2,1) = count;
    count = count + 1;
end

% Tower
% -----

Ntow = zeros(nNtow,1);                                                                                              % tower node ID array
Ntow(1,1) = Ntp(nNtop / 2 + 2,1);

for k = 1:nNtow - 1
    Nodes(count,1) = count;
    Nodes(count,2) = 0;
    Nodes(count,3) = 0;
    Nodes(count,4) = htot + k * htow / nNtow;
    Ntow(k + 1,1) = count;
    count = count + 1;
end

Ntt = max(Ntow);                                                                                                    % Tower top mass ID for locating RNA mass



% Foundation node 
% ---------------

Nodes(count,:) = [count, -b(1,1)/2-1, 0 , 0];   count = count + 1;
Nodes(count,:) = [count, -b(1,1)/2  , 0 , -1];  count = count + 1;
Nodes(count,:) = [count,  b(1,1)/2+1, 0 , 0];   count = count + 1;
Nodes(count,:) = [count,  b(1,1)/2  , 0 , -1];  count = count + 1;
Nfm = [Ntt+1:count]';                                                                                               % Foundation node ID array


% Last node == Reference node 
% --------------

Nodes(count,1) = count;
Nodes(count,2) = bbot + 1;
Nodes(count,3) = 0;
Nodes(count,4) = htot + 1;

% Node matrix
% ----------

Nodes = Nodes(1:count,:);

%% Element definition

% Element matrix
% --------------        
% Elements  = [EltID TypID SecID MatID n1 n2 n3]
% - EltID: elmement ID
% Beam element coordinate system: x-axis defined by n1 and n2,
%                                 xy-plane defined by n1, n2 and n3

Elements = zeros(500, 7);
Diameter = zeros(500, 3);
countEl = 1;                                                                                                        % counter representing element ID

% Botom member element definition
% -------------------------------

for k = 1:size(Nbm) - 1
    Elements(countEl,1) = countEl;
    Elements(countEl,2) = 1;
    Elements(countEl,3) = 3;
    Elements(countEl,4) = 1;
    Elements(countEl,5) = Nbm(k);       % n1
    Elements(countEl,6) = Nbm(k + 1);   % n2
    Elements(countEl,7) = count;        % n3
    
    Diameter(countEl,1) = countEl;
    Diameter(countEl,2) = Nbm(k);
    Diameter(countEl,3) = D(3,1);
    countEl = countEl + 1;
end


% Legs and braces element definition
% ----------------------------------

for j = 1:nseg
    
    % Leg element definition
    for k = 1:nNseg
        if j == 1
            Elements(countEl,1) = countEl;
            Elements(countEl,2) = 1;
            Elements(countEl,3) = 1;
            Elements(countEl,4) = 1;
            Elements(countEl,5) = Nleg1(k,j);
            Elements(countEl,6) = Nleg1(k + 1,j);
            Elements(countEl,7) = count;
            
            Diameter(countEl,1) = countEl;
            Diameter(countEl,2) = Nleg1(k,j);
            Diameter(countEl,3) = D(1,1);
            countEl = countEl + 1;
        else
            Elements(countEl,1) = countEl;
            Elements(countEl,2) = 1;
            Elements(countEl,3) = 2;
            Elements(countEl,4) = 1;
            Elements(countEl,5) = Nleg1(k,j);
            Elements(countEl,6) = Nleg1(k + 1,j);
            Elements(countEl,7) = count;
            
            Diameter(countEl,1) = countEl;
            Diameter(countEl,2) = Nleg1(k,j);
            Diameter(countEl,3) = D(2,1);
            countEl = countEl + 1;
        end
    end

    for k = 1:nNseg
        if j == 1
            Elements(countEl,1) = countEl;
            Elements(countEl,2) = 1;
            Elements(countEl,3) = 1;
            Elements(countEl,4) = 1;
            Elements(countEl,5) = Nleg2(k,j);
            Elements(countEl,6) = Nleg2(k + 1,j);
            Elements(countEl,7) = count;
            
            Diameter(countEl,1) = countEl;
            Diameter(countEl,2) = Nleg2(k,j);
            Diameter(countEl,3) = D(1,1);
            countEl = countEl + 1;
        else
            Elements(countEl,1) = countEl;
            Elements(countEl,2) = 1;
            Elements(countEl,3) = 2;
            Elements(countEl,4) = 1;
            Elements(countEl,5) = Nleg2(k,j);
            Elements(countEl,6) = Nleg2(k + 1,j);
            Elements(countEl,7) = count;
            
            Diameter(countEl,1) = countEl;
            Diameter(countEl,2) = Nleg2(k,j);
            Diameter(countEl,3) = D(2,1);
            countEl = countEl + 1;
        end
    end
    
    % Brace element definition
    for k = 1:nNseg
        Elements(countEl,1) = countEl;
        Elements(countEl,2) = 1;
        Elements(countEl,3) = 3;
        Elements(countEl,4) = 1;
        Elements(countEl,5) = Nbr1(k,j);
        Elements(countEl,6) = Nbr1(k + 1,j);
        Elements(countEl,7) = count;
        
        Diameter(countEl,1) = countEl;
        Diameter(countEl,2) = Nbr1(k,j);
        Diameter(countEl,3) = D(3,1);
        countEl = countEl + 1;
    end

    for k = 1:nNseg
        Elements(countEl,1) = countEl;
        Elements(countEl,2) = 1;
        Elements(countEl,3) = 3;
        Elements(countEl,4) = 1;
        Elements(countEl,5) = Nbr2(k,j);
        Elements(countEl,6) = Nbr2(k + 1,j);
        Elements(countEl,7) = count;
        
        Diameter(countEl,1) = countEl;
        Diameter(countEl,2) = Nbr2(k,j);
        Diameter(countEl,3) = D(3,1);
        countEl = countEl + 1;
    end
end

% Transition piece element definition
% ----------------

for k = 1:size(Ntp) - 1
    Elements(countEl,1) = countEl;
    Elements(countEl,2) = 1;
    Elements(countEl,3) = 4;
    Elements(countEl,4) = 2;
    Elements(countEl,5) = Ntp(k);
    Elements(countEl,6) = Ntp(k + 1);
    Elements(countEl,7) = count;
    
    Diameter(countEl,1) = countEl;
    Diameter(countEl,2) = Ntp(k);
    Diameter(countEl,3) = NaN;
    countEl = countEl + 1;
end

% Tower element definition
% ------------------------

    for k = 1:size(Ntow) - 1
        Elements(countEl,1) = countEl;
        Elements(countEl,2) = 1;
        Elements(countEl,3) = k + 4;
        Elements(countEl,4) = 1;
        Elements(countEl,5) = Ntow(k);
        Elements(countEl,6) = Ntow(k + 1);
        Elements(countEl,7) = count;
        
        Diameter(countEl,1) = countEl;
        Diameter(countEl,2) = Ntow(k);
        Diameter(countEl,3) = D(k + 3,1);
        countEl = countEl + 1;
    end
 
    
% Foundation member element definition
% -------------------------------

Elements(countEl,:)  = [countEl 2 sfnd(1) 3  36  1 count]; countEl = countEl + 1;
Elements(countEl,:)  = [countEl 2 sfnd(2) 3  37  1 count]; countEl = countEl + 1;
Elements(countEl,:)  = [countEl 2 sfnd(1) 3  38  2 count]; countEl = countEl + 1;
Elements(countEl,:)  = [countEl 2 sfnd(2) 3  39  2 count]; countEl = countEl + 1;    

% Element matrix
% ----------

Elements = Elements(1:countEl - 1,:);
Diameter = Diameter(1:countEl - 1,:);

%% title('Finite element model');

%% DOF definition
% DOF -> [2.01 2.02 2.06 3.01 ...].'
%   2.01 represents the 1st DOF (UX) of node 2
%   5.06 represents the 6st DOF (ROTZ) of node 5
%   0.04 represents the 4th DOF (ROTX) of all nodes
%   1.00 represents all DOFs of node 1

% Assemble DOF matrix
% -------------------
% Assemble a column matrix containing all DOFs at which stiffness is
% present in the model:

DOF=getdof(Elements,Types);

% Remove redundant DOFs
% ---------------------
% With foundation this changes to: 
% Remove all DOFs equal to zero from the vector:
%  - 2D analysis: select only UX,UZ,ROTY
%  - clamp node 36
%  - clamp node 37
%  - clamp node 38
%  - clamp node 39

seldof=[0.02; 0.04; 0.06; 36.00; 37.00; 38.00; 39.00];
DOF=removedof(DOF,seldof);

clear seldof

%% System matrix assembly

% Assembly of K and M
%     K: stiffness matrix
%     M: consistent mass matrix

[K,M]=asmkm(Nodes,Elements,Types,Sections,Materials,DOF);

% Add top mass
tmDOF = Ntt + 0.01;                     % [Ntt + 0.01, Ntt + 0.03;]
[~,Mtm] = selectdof(DOF,tmDOF);
M(Mtm,Mtm) = M(Mtm,Mtm) + Mtop;

%% Eigensolution

% Solution of the eigenvalue problem K * phi = Lambda * M * phi
% --------------
% Returns the eigenfrequencies in [rad/s] (column matrix Omega) and the
% mass-normalized mode shapes (Phi, every column contains 1 mode shape).

[Phi,Omega]=eigfem(K,M);

%% Visualisations
% % comment this section if function is invoked multiple times:
% 
% % Node and element plot
% % --------------
% 
% figure(); hold on;
% plotnodes(Nodes);
% plotelem(Nodes,Elements,Types,'Numbering','off');
% 
% % Eigenmode plot
% % Animate the N eigenmodes: 
% % --------------
% 
% for ind = 1:5
%     figure;
%     animdisp(Nodes,Elements,Types,DOF,Phi(:,ind));
%     title(sprintf('Mode %d',ind));
% end
% 
% % Natural frequencies
% % -------------------
% 
% disp('Eigenfrequencies [Hz]:');
% disp(Omega(1:20)/2/pi);