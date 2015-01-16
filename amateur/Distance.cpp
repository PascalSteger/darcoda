#include "Global.h"
#include "Distance.h"

Distance::Distance( real r, unsigned i ){
  r_ = r;
  i_ = i;
}

Distance::Distance( const Distance& d ){
  r_ = d.getr();
  i_ = d.geti();
}

bool
Distance::operator<( const Distance& d ) const {
  return r_ < d.getr();
}

real
Distance::getr( void ) const {
  return r_;
}

unsigned
Distance::geti( void ) const {
  return i_;
}
