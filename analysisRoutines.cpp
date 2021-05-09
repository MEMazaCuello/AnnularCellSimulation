/**
  * "analysisRoutines.cpp":
  * --------------------
  * Implementation of methods for analysis of the configurations.
  * For declaration details, see "analysisRoutines.h".
  *
  * Needs: "AnnularCell.h".
  *
  * --------------------
  * Includes methods for studying the
  *   a) Order parameters,
  *   b) Clusters, and
  *   c) Defects of tetratic field.
  *
  * --------------------
  * Last modified: 2021-05-09
  * By: M. E. Maza-Cuello
  */

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "analysisRoutines.h"

// External parameters from "parameters.cpp"
extern const double PI;
extern const double HALF_PI;
extern const double WIDTH;
extern const double LENGTH;
extern const double DIAGONAL;
extern const double AVG_RADIUS;
extern const int    DEFECT_MIN_SIZE;

// Order Parameters analysis __________________________________________

std::vector<int> getAveragingIndexes(const AnnularCell& cell, const Rod& rod)
{
  std::vector<int> neighbors;

  for (int i = 0; i < NUMBER_OF_RODS; i++)
  {
    if ( cell.getRod(i).isWithinRadius(rod,AVG_RADIUS) )
    {
        neighbors.push_back(i);
    }
  }

  return neighbors;
}

std::vector<int> getAveragingIndexes(const AnnularCell& cell, const int& index)
{
  return getAveragingIndexes(cell, cell.getRod(index));
}

std::vector<double> getQ4(AnnularCell& cell)
{
  double tilt;
  double zeta;
  double sum_cos4;
  double sum_sin4;

  std::vector<int> indexes;
  std::vector<double> q4(cell.m_NUMBER_OF_PARTICLES);

  for (int i = 0; i < NUMBER_OF_RODS; i++)
  {
    // Neighboring indexes to perform local average
    indexes = getAveragingIndexes(cell, i);

    // Getting the local nematic eigen direction "tilt"
    sum_cos4 = 0.0d;
    sum_sin4 = 0.0d;
    for (int index : indexes)
    {
      zeta = 2.0d*cell.m_bundle[index].m_angle;
      sum_cos4 += std::cos(zeta);
      sum_sin4 += std::sin(zeta);
    }
    tilt = 0.5d*std::atan2(sum_sin4,sum_cos4);

    // Compute local tetratic order parameter (q4) relative to "tilt"
    sum_cos4 = 0.0d;
    sum_sin4 = 0.0d;
    for (int index : indexes)
    {
      zeta = 4.0d*(cell.m_bundle[index].m_angle-tilt);
      sum_cos4 += std::cos(zeta);
      sum_sin4 += std::sin(zeta);
    }

    q4[i] = std::sqrt(sum_cos4*sum_cos4 + sum_sin4*sum_sin4)/indexes.size();
  }

  return q4;
}

std::vector<double> getQ4FromFile(const std::string& filepath, const int& numRodsInFile)
{
  double dummy;
  std::ifstream savedconfiguration;
  savedconfiguration.open(filepath);

  std::vector<double> q4(numRodsInFile);

  for (int i = 0; i < numRodsInFile; i++)  // From order parameters files
  {
    savedconfiguration >> dummy >> dummy >> dummy >> dummy
                       >> dummy >> dummy >> q4[i] >> dummy;
  }
  savedconfiguration.close();

  return q4;
}

void getOrderParameters(AnnularCell& cell, std::string& filename)
{
  static double scale = 2.0d*PI/(1.2d*LENGTH);
  static std::ofstream Qdat, Qcsv;

  Qdat.open(filename+".txt");
  Qcsv.open(filename+".csv");
  Qcsv << "index,x,y,angle,tilt,q2,q4,qs" << std::endl;

  double tilt, zeta;
  double sum_cos2, sum_sin2;
  double sum_cos4, sum_sin4;
  double sum_cosS, sum_sinS;

  std::vector<int> indexes;
  for (int i = 0; i < NUMBER_OF_RODS; i++)
  {
    sum_cos2 = 0.0d;
    sum_sin2 = 0.0d;

    indexes = getAveragingIndexes(cell, i);

    // Getting the local nematic eigen direction "tilt"
    for (int index : indexes)
    {
      sum_cos2 += std::cos(2.0d*cell.m_bundle[index].m_angle);
      sum_sin2 += std::sin(2.0d*cell.m_bundle[index].m_angle);
    }
    tilt = 0.5d*std::atan2(sum_sin2,sum_cos2);

    // Compute local order parameters (q2, q4, qs) relative to "tilt"
    sum_cos2 = 0.0d; sum_sin2 = 0.0d;
    sum_cos4 = 0.0d; sum_sin4 = 0.0d;
    sum_cosS = 0.0d; sum_sinS = 0.0d;

    for (int index : indexes)
    {
      // Local nematic order parameter q2
      zeta = 2.0d*(cell.m_bundle[index].m_angle-tilt);
      sum_cos2 += std::cos(zeta);
      sum_sin2 += std::sin(zeta);

      // Local tetratic order parameter q4
      sum_cos4 += std::cos(2.0d*zeta);
      sum_sin4 += std::sin(2.0d*zeta);

      // Local smectic order parameter qs
      zeta = scale*(  std::cos(tilt)*(cell.m_bundle[index].m_xPos-cell.m_bundle[i].m_xPos)
                    + std::sin(tilt)*(cell.m_bundle[index].m_yPos-cell.m_bundle[i].m_yPos) );
      sum_cosS += std::cos(zeta); // cos( (2pi/(1.2*LENGTH)) * zeta)
      sum_sinS += std::sin(zeta);
    }

    zeta = 1.0d/indexes.size();
    Qdat << i+1 << " " << cell.m_bundle[i].m_xPos << " " << cell.m_bundle[i].m_yPos << " "
         << cell.m_bundle[i].m_angle << " " << tilt << " "                           // Save configuration (and tilt (eigen)angle)
         << std::sqrt(sum_cos2*sum_cos2 + sum_sin2*sum_sin2)*zeta << " "             // Nematic order parameter Q1
         << std::sqrt(sum_cos4*sum_cos4 + sum_sin4*sum_sin4)*zeta << " "             // Tetratic order parameter Q2
         << std::sqrt(sum_cosS*sum_cosS + sum_sinS*sum_sinS)*zeta << std::endl;      // Smectic order parameter Qs

    Qcsv << i+1 << "," << cell.m_bundle[i].m_xPos << "," << cell.m_bundle[i].m_yPos << ","
         << cell.m_bundle[i].m_angle << "," << tilt << ","                           // Save configuration (and tilt (eigen)angle)
         << std::sqrt(sum_cos2*sum_cos2 + sum_sin2*sum_sin2)*zeta << ","             // Nematic order parameter Q1
         << std::sqrt(sum_cos4*sum_cos4 + sum_sin4*sum_sin4)*zeta << ","             // Tetratic order parameter Q2
         << std::sqrt(sum_cosS*sum_cosS + sum_sinS*sum_sinS)*zeta << std::endl;      // Smectic order parameter Qs
  }
  Qdat.close();
  Qcsv.close();
}

// Auxiliary tree-structure ___________________________________________

std::map< int, std::vector<int> > getTrees(std::vector<Link> links)
{
  // Links must be partially ascending-ordered w.r.t link.left
  // Hence "trees" will be a map with ordered keys.
  // Each key is the root of a tree-structure representing one tree.
  std::map< int, std::vector<int> > trees;

  // Roots are saved as non-positive (<= 0) integers
  int root;
  int leaf;

  bool leftHasRoot;
  bool rightHasRoot;

  for(Link lk : links)
  {
    leftHasRoot  = trees.count(lk.left);
    rightHasRoot = trees.count(lk.right);

    if ( !rightHasRoot && !leftHasRoot )
    {
      // New cluster
      trees[lk.left].push_back(lk.left);

      trees[lk.left].push_back(lk.right);  // Root: lk.left
      trees[lk.right].push_back(-lk.left); // Leaf: lk.right
    }
    else if ( !rightHasRoot && leftHasRoot )
    {
      if ( trees[lk.left][0] < lk.left ) // lk.left is not its own root
      {
        root = -trees[lk.left][0];
        trees[root].push_back(lk.right);   // Root: root
        trees[lk.right].push_back(-root);  // Leaf: lk.right
      }
      else // lk.left is its own root
      {
        trees[lk.left].push_back(lk.right);  // Root: lk.left
        trees[lk.right].push_back(-lk.left); // Leaf: lk.right
      }
    }
    else if ( rightHasRoot && !leftHasRoot )
    {
      // "root" will be non-negative (>= 0) because Links are ordered
      root = -trees[lk.right][0];

      // root < lk.left because Links are ordered
      trees[root].push_back(lk.left);  // Root: root
      trees[lk.left].push_back(-root); // Leaf: lk.left
    }
    else // rightHasRoot && leftHasRoot
    {
      // Merging of trees
      if ( trees[lk.left][0] < lk.left ) // lk.left is not its own root
      {
        if ( trees[lk.left][0] < trees[lk.right][0] ) // left into right
        {
          // "root" will be non-negative (>= 0) because Links are ordered
          root = -trees[lk.right][0];
          leaf = -trees[lk.left][0];
          for (int i : trees[leaf])
          {
            trees[i][0] = (-root);  // Root: root
            trees[root].push_back(i);
          }
          trees[leaf].resize(1);
        }
        else if ( trees[lk.right][0] < trees[lk.left][0] ) // right into left
        {
          // "root" will be non-negative (>= 0) because Links are ordered
          root = -trees[lk.left][0];
          leaf = -trees[lk.right][0];
          for (int i : trees[leaf])
          {
            trees[i][0] = (-root);  // Root: root
            trees[root].push_back(i);
          }
          trees[leaf].resize(1);
        }
      }
      else // lk.left is its own root
      {
        // "root" will be non-negative (>= 0) because Links are ordered
        root = -trees[lk.right][0];
        leaf = lk.left;
        for (int i : trees[leaf])
        {
          trees[i][0] = (-root);  // Root: root
          trees[root].push_back(i);
        }
        trees[leaf].resize(1);
      }
    }
  }

  return trees;
}

void eraseLeaves(std::map< int, std::vector<int>>& trees, const int& minSize)
{
  for (auto const& tree : trees)
  {
    if (tree.second.size() < minSize)
    {
      trees.erase(tree.first);
    }
  }
}

// Clusters analysis __________________________________________________

bool areInSameCluster(Rod self, Rod other)
{
  static double maxAngle = PI/18.0d; // 10 degrees
  static double maxDist  = 1.8d*WIDTH;

  // Distance criterion
  if  ( !self.isWithinRadius(other,maxDist) )
  {
    return false;
  }

  // Orientation criterion
  double angleDiff = std::abs(self.m_angle-other.m_angle);
  if (angleDiff > HALF_PI)
  {
    angleDiff = PI - angleDiff;
  }

  return (angleDiff < maxAngle);
}

std::vector<Link> getClusterLinks(AnnularCell& cell)
{
  std::vector<Link> links;
  Link link;

  // Convention: link.left < link.right
  // Thus "links" is a partially ordered vector w.r.t. link.left.
  for (int i = 0; i < cell.m_NUMBER_OF_PARTICLES; i++)
  {
    cell.m_grid.setNeighbors(cell.m_grid.getCoords(cell.m_bundle[i]));
    for (int index : cell.m_grid.m_neighbors)
    {
      if (i < index && areInSameCluster(cell.getRod(index), cell.getRod(i)))
      {
        link.left  = i;
        link.right = index;
        links.push_back(link);
      }
    }
  }

  return links;
}

std::map< int, std::vector<int> > getClusters(AnnularCell& cell)
{
  return getTrees(getClusterLinks(cell));
}

std::vector<int> getClusterLabels(AnnularCell& cell)
{
  std::map< int, std::vector<int> > clusters = getClusters(cell);

  // "label[index] = -1" means that particle "index" is isolated.
  std::vector<int> labels(cell.m_NUMBER_OF_PARTICLES,-1);

  for (auto const& cluster : clusters)
  {
    if (cluster.second[0] >= 0)
    {
      for (int const& i : cluster.second)
      {
        labels[i] = cluster.first;
      }
    }
  }

  return labels;
}

void saveCellWithClusters(AnnularCell& cell, std::string& filename)
{
  static std::ofstream Qdat, Qcsv;

  Qdat.open(filename+".txt");
  Qcsv.open(filename+".csv");
  Qcsv << "index,x,y,angle,cluster" << std::endl;

  std::vector<int> clusters = getClusterLabels(cell);

  for(int i = 0; i < NUMBER_OF_RODS; i++)
  {
    Qdat << i+1 << " " << cell.m_bundle[i].m_xPos << " " << cell.m_bundle[i].m_yPos << " "
         << cell.m_bundle[i].m_angle << " " << clusters[i] << std::endl;
    Qcsv << i+1 << "," << cell.m_bundle[i].m_xPos << "," << cell.m_bundle[i].m_yPos << ","
         << cell.m_bundle[i].m_angle << "," << clusters[i] << std::endl;
  }
  Qdat.close();
  Qcsv.close();
}

// Defects analysis ___________________________________________________

std::vector<int> getParticlesInDefects(const std::vector<double>& q4)
{
  // Threshold value for considering a particle part of a defect
  static const double minQ4value = 0.4;

  std::vector<int> particles;

  for (int i = 0; i < q4.size() ; i++)
  {
      if (q4[i] < minQ4value)
      {
        particles.push_back(i);
      }
  }

  return particles;
}

bool areInSameDefect(Rod self, Rod other)
{
  static const double maxDist = 2.0d*DIAGONAL;

  return self.isWithinRadius(other,maxDist);
}

std::vector<Link> getDefectsLinks(AnnularCell& cell, const std::vector<double>& q4)
{
  std::vector<Link> links;
  Link link;

  // Particles that belong to defects
  std::vector<int> particles = getParticlesInDefects(q4);

  // Convention: link.left < link.right
  // Thus "links" is a partially ordered vector w.r.t. link.left.
  for (int i : particles)
  {
    for (int index : particles)
    {
      if (i < index && areInSameDefect(cell.getRod(index), cell.getRod(i)))
      {
        link.left  = i;
        link.right = index;
        links.push_back(link);
      }
    }
  }

  return links;
}

inline std::vector<Link> getDefectsLinks(AnnularCell& cell)
{
  return getDefectsLinks(cell,getQ4(cell));
}

std::map< int, std::vector<int> > getDefects(AnnularCell& cell, const std::vector<double>& q4)
{
  return getTrees(getDefectsLinks(cell,q4));
}

inline std::map< int, std::vector<int> > getDefects(AnnularCell& cell)
{
  return getTrees(getDefectsLinks(cell));
}

std::vector<int> getDefectsLabels(AnnularCell& cell, const std::vector<double>& q4)
{
  std::map< int, std::vector<int> > defects = getDefects(cell,q4);

  // "label[index] = -1" means that particle "index" is isolated.
  std::vector<int> labels(cell.m_NUMBER_OF_PARTICLES,-1);

  for (auto const& defect : defects)
  {
    if (defect.second[0] >= 0)
    {
      for (int const& i : defect.second)
      {
        labels[i] = defect.first;
      }
    }
  }

  return labels;
}

inline std::vector<int> getDefectsLabels(AnnularCell& cell)
{
  return getDefectsLabels(cell,getQ4(cell));
}

std::map<int,int> getDefectsSizes(std::map< int, std::vector<int> > defects)
{
  eraseLeaves(defects,DEFECT_MIN_SIZE);

  std::map<int,int> sizes;
  for (auto const& defect: defects)
  {
    sizes[defect.first] = defect.second.size();
  }

  return sizes;
}

std::map<int,std::vector<double>> getDefectsCartesianCoords(std::map<int,std::vector<int>> defects, AnnularCell& cell)
{
  eraseLeaves(defects,DEFECT_MIN_SIZE);

  std::map<int,std::vector<double>> positions;

  double x;
  double y;
  Rod rod;
  for (auto const& defect: defects)
  {
    x = 0;
    y = 0;
    for (int i : defect.second)
    {
      rod = cell.getRod(i);
      x += rod.m_xPos;
      y += rod.m_yPos;
    }
    x = x/defect.second.size();
    y = y/defect.second.size();

    positions[defect.first] = std::vector<double> {x,y};
  }

  return positions;
}

std::map<int,std::vector<double>> getDefectsPolarCoords(std::map<int,std::vector<int>> defects, AnnularCell& cell)
{
  std::map<int,std::vector<double>> positions = getDefectsCartesianCoords(defects,cell);

  double x;
  double y;
  for (auto & position : positions)
  {
    x = position.second[0];
    y = position.second[1];

    position.second[0] = std::sqrt(x*x + y*y);
    position.second[1] = std::atan2(y,x);
  }

  return positions;
}

void saveCellWithDefects(AnnularCell& cell, std::string& filename, const std::vector<double>& q4)
{
  static std::ofstream Qdat, Qcsv;

  Qdat.open(filename+".txt");
  Qcsv.open(filename+".csv");
  Qcsv << "index,x,y,angle,defect" << std::endl;

  std::vector<int> defects = getDefectsLabels(cell,q4);

  for (int i = 0; i < NUMBER_OF_RODS; i++)
  {
    Qdat << i+1 << " " << cell.m_bundle[i].m_xPos << " " << cell.m_bundle[i].m_yPos << " "
         << cell.m_bundle[i].m_angle << " " << defects[i] << std::endl;
    Qcsv << i+1 << "," << cell.m_bundle[i].m_xPos << "," << cell.m_bundle[i].m_yPos << ","
         << cell.m_bundle[i].m_angle << "," << defects[i] << std::endl;
  }
  Qdat.close();
  Qcsv.close();
}

inline void saveCellWithDefects(AnnularCell& cell, std::string& filename)
{
  saveCellWithDefects(cell, filename, getQ4(cell));
}

std::vector<int> getHistogramCounts(std::vector<double> v, double inf, double sup, int bins)
{
  std::vector<int> hist(bins,0);
  double delta = (sup-inf)/bins;

  int bin;
  for (double const& val : v)
  {
    bin = int(floor((val - inf)/delta));
    try
    {
      hist.at(bin)++;
    }
    catch (const std::exception& e)
    {
      std::cout << "Value: " << val <<
                   " is outside the interval [ " << inf << " , " << sup << " ]"
                << " and thus not included in the histogram." << std::endl;
    }
  }

  return hist;
}
