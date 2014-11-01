#ifndef DISTANCE_H
#define DISTANCE_H

#include "Global.h"

class Distance {
 private:
  real r_;
  unsigned i_;

 public:
  Distance( real r, unsigned i );
  Distance( const Distance& d );
  bool operator<( const Distance& d ) const;
  real getr( void ) const;
  unsigned geti( void ) const;
};

#endif
