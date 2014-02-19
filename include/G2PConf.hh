// -*- C++ -*-

/* Config class
 * Use libconfig, a 3rd party package, to do the parsing.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Feb 2014, C. Gu, Modified for G2PRec.
//

#ifndef G2P_CONF_H
#define G2P_CONF_H

#include "libconfig.h++"

using namespace libconfig;

extern class Config* gConfig;

#endif