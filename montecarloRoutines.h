#ifndef MONTECARLOROUTINES_H_INCLUDED
#define MONTECARLOROUTINES_H_INCLUDED

#include "AnnularCell.h"

void stepMontecarlo(AnnularCell& cell);
void getOrderParameters(AnnularCell& cell);
void getOrderParameters(AnnularCell& cell, std::string& filename);

#endif // MONTECARLOROUTINES_H_INCLUDED
