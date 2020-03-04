// Gmsh project created on Thu Dec 06 01:42:08 2018
SetFactory("OpenCASCADE");

// Declaring mesh parameters
hp = 0.02;
L = 2;
ns = 5;

//Declare base points and first line (front of plate)
Point(2) = {0, hp/2, 0, 1.0};
Point(3) = {0, -hp/2, 0, 1.0};
Line(1) = {3, 2};

// Extrude from from of plate to create the plate
Extrude {L, 0, 0} {
  Line{1}; 
}//+
Transfinite Line {3, 2} = ns*50 Using Progression 0.995;
//+
Transfinite Line {1, 4} = ns Using Progression 1;
//+
Transfinite Surface {1};
//+
Physical Line("plate_free") = {3, 4, 2};
//+
Physical Line("plate_fixed") = {1};
//+
Physical Surface("plate") = {1};
