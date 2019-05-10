#include <cmath>

/** PRIMARY PARAMETERS **/
extern const double WIDTH    = 1.0d;
extern const double LENGTH   = 4.0d*WIDTH;  // LENGTH GREATER (OR EQUAL) THAN WIDTH!!
                                            // UNKNOWN BEHAVIOUR OTHERWISE
extern const double DIAGONAL = std::sqrt(WIDTH*WIDTH + LENGTH*LENGTH);

extern const double OUTER_RADIUS = 75.0d*WIDTH;
extern const double INNER_RADIUS = 10.0d*WIDTH; // INNER SMALLER THAN OUTER!!
                                                // UNKNOWN BEHAVIOUR OTHERWISE

extern const int    NUMBER_OF_RODS = 3255; // Desired packing fraction must be computed beforehand

/** AUXILIARY PARAMETERS **/
extern const double HALF_WIDTH    = 0.5d*WIDTH;
extern const double HALF_LENGTH   = 0.5d*LENGTH;
extern const double HALF_DIAGONAL = 0.5d*DIAGONAL;
extern const double ALPHA = std::atan2(WIDTH,LENGTH);  // Inner angle of the rectangles

extern const double HALF_PI    = std::asin(1.0d);
extern const double PI         = std::acos(-1.0d);
extern const double QUARTER_PI = 0.5d*HALF_PI;

extern const double PACKING_FRACTION = NUMBER_OF_RODS*WIDTH*LENGTH/(PI*(OUTER_RADIUS*OUTER_RADIUS-INNER_RADIUS*INNER_RADIUS));
