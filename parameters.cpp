/**
  * "parameters.cpp":
  * --------------------
  * Sets global simulation parameters used by other files.
  *
  * --------------------
  * Primary parameters:
  * Set characteristic lengths of the rods and the annular cell radii,
  * and the number of rods to be included in the cell.
  * They also include parameters relating to the Montecarlo simulation and
  * the analysis of clusters and defects within the cell.
  * Note: all lengths are defined with respect to the "WIDTH" of the rods.
  * Note: it is everywhere assumed that
  *   a) WIDTH <= LENGTH,
  *   b) INNER_RADIUS < OUTER_RADIUS,
  *   c) the specified NUMBER_OF_RODS can indeed fit inside the annular cell
  *      with given dimensions.
  * It is recommended to work with packing fractions lesser than 0.75.
  *
  * Auxiliary parameters:
  * Derived from primary parameters, used for optimization purposes.
  * Note: these parameters do not require user input.
  *
  * --------------------
  * Last modified: 2021-05-09
  * By: M. E. Maza-Cuello
  */

#include <cmath>

// Primary parameters _________________________________________________

// Rod
extern const double WIDTH  = 1.0d;
extern const double LENGTH = 4.0d*WIDTH;

// Annular Cell
extern const double OUTER_RADIUS = 70.0d*WIDTH;
extern const double INNER_RADIUS = 20.0d*WIDTH;
extern const int    NUMBER_OF_RODS = 2651;

// Montecarlo parameters

extern const int    SEED = 653248971;

/**
  * Set desired acceptance probability (in %) of a Montecarlo step.
  * A Montecarlo step consists of trying to slightly move each rod once.
  */
extern const double ACCEPTANCE = 0.25d;

// Order parameters

/**
  * Set desired AVG_RADIUS of circular region around rod to obtain
  * local averages of order parameters.
  */
extern const double AVG_RADIUS = 4.0d*LENGTH;

// Clusters

/**
  * Set minimum number of rods to be considered a cluster
  */
extern const int CLUSTER_MIN_SIZE = 2;

// Defects

/**
  * Set minimum number of rods to be considered a defect
  */
extern const int DEFECT_MIN_SIZE = 5;

// Auxiliary parameters _______________________________________________

// Rod
extern const double DIAGONAL = std::sqrt(WIDTH*WIDTH + LENGTH*LENGTH);
extern const double HALF_WIDTH    = 0.5d*WIDTH;
extern const double HALF_LENGTH   = 0.5d*LENGTH;
extern const double HALF_DIAGONAL = 0.5d*DIAGONAL;
extern const double ALPHA = std::atan2(WIDTH,LENGTH);  // Inner angle of the rectangles

// Angles
extern const double HALF_PI = std::asin(1.0d);
extern const double PI      = std::acos(-1.0d);

// Annular Cell
extern const double PACKING_FRACTION = NUMBER_OF_RODS*WIDTH*LENGTH/(PI*(OUTER_RADIUS*OUTER_RADIUS-INNER_RADIUS*INNER_RADIUS));
