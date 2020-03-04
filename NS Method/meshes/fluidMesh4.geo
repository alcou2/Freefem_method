// Gmsh for simple flag problem

SetFactory("OpenCASCADE");


// Declaring mesh parameters
hp = 0.02;
L = 2;
h = 0.5;
RL1 = 0.5;
RL21 = 0.5;
RL22 = 2.5;
ns = 5;
nf = ns*2;



//Declare base points and first line (front of plate)
Point(2) = {0, hp/2, 0, 1.0};
Point(3) = {0, -hp/2, 0, 1.0};
Line(1) = {3, 2};

// Extrude from front of plate to create the plate
Extrude {L, 0, 0} {
  Line{1}; 
}

// Extrude to create fluid domain above plate
Extrude {0, (h-hp)*0.7, 0} {
  Line{3}; 
}

// Extrude to create fluid domain under plate
Extrude {0, -(h-hp)*0.7, 0} {
  Line{2}; 
}

// Extrude to create fluid domain in front of the plate
Extrude {-L*RL1, 0, 0} {
  Line{1}; Line{5}; Line{8}; 
}

// Extrude to create fluid domain behind the plate (close)
Extrude {L*RL21, 0, 0} {
  Line{6}; Line{4}; Line{9}; 
}

// Extrude to create fluid domain behind the plate (far)
Extrude {L*RL22, 0, 0} {
  Line{20}; Line{22}; Line{24}; 
}

// Extrude to create fluid domain above plate (close to top)
Extrude {0, (h-hp)*0.3, 0} {
  Line{14}; Line{7}; Line{19}; Line{26};
}

// Extrude to create fluid domain bellow the plate (close to bottom)
Extrude {0, -(h-hp)*0.3, 0} {
  Line{16}; Line{10}; Line{23}; Line{30};
}

// Discretize the diffent lines


// All vectical lines above plate (far)
Transfinite Line {33, 32, 35, 37, 39} = nf Using Progression 0.85;

// All vectical lines above plate (close)
Transfinite Line {15, 5, 6, 20, 27} = nf*2 Using Progression 1.15;

// All vertical lines at the plate height
Transfinite Line {4, 22, 29} = ns Using Progression 1;
//Transfinite Line {22, 29} = 3 Using Progression 1;

// All vectical lines below plate (close)
Transfinite Line {17, 8, 9, 24, 31} = nf*2 Using Progression 1.15;

// All vectical lines below plate (far)
Transfinite Line {42, 41, 44, 46, 48} = nf Using Progression 0.85;



// Horizontal lines in front of plate
Transfinite Line {43, 16, 11, 12, 14, 34} = nf*4 Using Progression 1.03;

// Horizontal lines above and below plate
Transfinite Line {45, 10, 2, 3, 7, 36} = 50*ns Using Progression 0.995;

// Horizontal lines fluid domain behind the plate (close)
Transfinite Line {47, 23, 21, 18, 19, 38} = nf*6 Using Progression 1.04;

// Horizontal lines fluid domain behind the plate (far)
Transfinite Line {49, 30, 28, 25, 26, 40} = nf*12 Using Progression 1;



// Declaring surfaces for structured mesh
Transfinite Surface {5} = {12, 6, 2, 11};
Transfinite Surface {2} = {6, 7, 5, 2};
Transfinite Surface {7} = {7, 15, 14, 5};
Transfinite Surface {10} = {15, 19, 18, 14};
Transfinite Surface {11} = {14, 18, 20, 16};
Transfinite Surface {12} = {16, 20, 21, 17};
Transfinite Surface {9} = {4, 16, 17, 9};
Transfinite Surface {3} = {3, 4, 9, 8};
Transfinite Surface {6} = {10, 3, 8, 13};
//Transfinite Surface {4} = {11, 2, 3, 10}; // Front of plate
//Transfinite Surface {1} = {2, 5, 4, 3};
Transfinite Surface {8} = {5, 14, 16, 4};

Transfinite Surface {13} = {23,22,6,12};
Transfinite Surface {14} = {22,24,7,6};
Transfinite Surface {15} = {24,25,15,7};
Transfinite Surface {16} = {25,26,19,15};

Transfinite Surface {17} = {13,8,27,28};
Transfinite Surface {18} = {8,9,29,27};
Transfinite Surface {19} = {9,17,30,29};
Transfinite Surface {20} = {17,21,31,30};


// Declare physical groups
Physical Line("plate_free") = {2, 4, 3};
Physical Line("plate_fixed") = {1};
Physical Line("inlet_up") = {15,33};
Physical Line("inlet_down") = {17,42};
Physical Line("outlet") = {48,31, 29, 27,39};
Physical Line("wall") = {34,36,38,40,43,45,47,49};
Physical Line("plate_inf") = {11,12};

Physical Surface("fluid") = {5, 2, 7, 10, 11, 12, 9, 8, 3,6,13,14,15,16,17,18,19,20};
//Physical Surface("plate") = {1};

//Recombine Surface {5, 2, 7, 10, 4, 1, 8, 11, 6, 3, 9, 12};

//Delete the plate from this fluid mesh
Recursive Delete {
  Surface{1};
  Surface{4}; 
}
