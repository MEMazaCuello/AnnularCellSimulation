/**** 
  'parameters.cpp':
  --------------------
  Sets global simulation parameters and lengths, called from other files.
  
  --------------------
  PRIMARY PARAMETERS define characteristic lengths of the rods,
  'WIDTH' and 'LENGTH'; annular cell radii, 'OUTER_RADIUS' and
  'INNER_RADIUS'; and 'NUMBER_OF_RODS' to be included in the cell.
  
  Note: all lengths are defined with respect to the 'WIDTH'.
  
  Note: the program expects LENGTH >= WIDTH 
                            INNER_RADIUS <= OUTER_RADIUS,
        and it is assumed that the specified NUMBER_OF_RODS can indeed 
        fit inside the annular cell with given dimensions.
  
  AUXILIARY PARAMETERS are computed from PRIMARY PARAMETERS values
  and are used to optimise certain parts of the code, in particular 
  for the detection of overlapping of rods between themselfes and with boundaries. 
  
  Note: these parameters are derived from previous ones and they 
        do not require user input.
        
  --------------------
  Last modified: 2019-05-12 
  By: M. E. Maza Cuello
****/

#include <cmath>

/** PRIMARY PARAMETERS **/
// Rods
extern const double WIDTH    = 1.0d;            // Set desired WIDTH of rods, used as unity of length
extern const double LENGTH   = 4.0d*WIDTH;      // Set desired LENGTH of rods
// Annular Cell
extern const double OUTER_RADIUS = 75.0d*WIDTH; // Set desired OUTER_RADIUS of the annular cell
extern const double INNER_RADIUS = 10.0d*WIDTH; // Set desired INNER_RADIUS of the annular cell
// Number of Rods
extern const int    NUMBER_OF_RODS = 3255;      // Set desired NUMBER_OF_RODS

/** AUXILIARY PARAMETERS **/
// Rods
extern const double DIAGONAL = std::sqrt(WIDTH*WIDTH + LENGTH*LENGTH);
extern const double ALPHA = std::atan2(WIDTH,LENGTH);  // Inner angle of the rectangles
extern const double HALF_WIDTH    = 0.5d*WIDTH;
extern const double HALF_LENGTH   = 0.5d*LENGTH;
extern const double HALF_DIAGONAL = 0.5d*DIAGONAL;
// Usefull angles (rads)
extern const double HALF_PI    = std::asin(1.0d);
extern const double PI         = std::acos(-1.0d);
extern const double QUARTER_PI = 0.5d*HALF_PI;
// Packing fraction for given parameters: area covered by rods over area of annular cell
extern const double PACKING_FRACTION = NUMBER_OF_RODS*WIDTH*LENGTH/(PI*(OUTER_RADIUS*OUTER_RADIUS-INNER_RADIUS*INNER_RADIUS));

