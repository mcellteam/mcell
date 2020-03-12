// some definitions shared among all MCell's components
// out of any namespace
// the reason to avoid namespaces is the capability to detach the BNG library
// from MCell
//
// TODO: use some master namespace? how should we call it?

#ifndef __SHARED_DEFINES_H__
#define __SHARED_DEFINES_H__

#include <stdint.h>
#include <climits>
#include <cfloat>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <map>
#include <unordered_map>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_set.hpp>

// ---------------------------------- float types ----------------------------------

#define FLOAT_T_BYTES 8

// float_t is also defined in mathdef.h, we need to enclose it into
// a namespace and then copy the definition in all define headers
namespace Common {

#if FLOAT_T_BYTES == 8
typedef double float_t; // will be changed to float
#else
#error "Base type float32 is not supported yet"
#endif

}


#if FLOAT_T_BYTES == 8
const Common::float_t EPS = 1e-12; // same as EPS_C
const Common::float_t SQRT_EPS = 1e-6;
const Common::float_t GIGANTIC4 = 1e140;
#else
#error "Base type float32 is not supported yet"
#endif

// ---------------------------------- fixed constants and specific typedefs -------------------
const Common::float_t POS_INVALID = FLT_MAX; // cannot be NAN because we cannot do any comparison with NANs

const Common::float_t TIME_INVALID = -1;
const Common::float_t TIME_FOREVER = FLT_MAX; // this max is sufficient for both float and double
const Common::float_t TIME_SIMULATION_START = 0;


const uint ID_INVALID = UINT32_MAX; // general invalid index, should not be used when a definition for a specific type is available
const uint INDEX_INVALID = UINT32_MAX; // general invalid index, should not be used when a definition for a specific type is available


// molecule id is a unique identifier of a molecule,
// no 2 molecules may have the same ID in the course of a simulation (at least for now)
// TODO: this should be something else, not ID, but neither index,
// it changes
typedef uint molecule_id_t;
const molecule_id_t MOLECULE_ID_INVALID = ID_INVALID;


#ifndef NDEBUG
// TODO: probably make this enabled only for Eclipse, we want the debug build to behave exactly as the release build
#define INDEXER_WA // Don't know yet how to convince Eclipse to correctly index boost containers
#endif

#if defined(NDEBUG) && defined(INDEXER_WA)
#warning "INDEXER_WA is enabled and this will lead to lower performance"
#endif


#ifndef INDEXER_WA
template<class T, class Allocator=boost::container::new_allocator<T>>
  using small_vector = boost::container::small_vector<T, 8, Allocator>;

template<class T, typename Compare = std::less<T>, class Allocator=boost::container::new_allocator<T>>
  using base_flat_set = boost::container::flat_set<T, Compare, Allocator>;
#else
template<typename T, typename _Alloc = std::allocator<T>  >
  using small_vector = std::vector<T, _Alloc>;

template<typename T, typename _Compare = std::less<T>, typename _Alloc = std::allocator<T>  >
  using base_flat_set = std::set<T, _Compare, _Alloc>;
#endif


/**
 * Template class used to hold sets of ids or indices of molecules or other items,
 * extended to check for unique insertions and checked removals.
 */
template<typename T>
class uint_set: public base_flat_set<T> {
public:
  // insert with check that the item is not there yet
  // for insertions without this check use 'insert'
  void insert_unique(const T id_or_index) {
    assert(this->count(id_or_index) == 0);
    this->insert(id_or_index);
  }

  // erase with check that the item is present
  // for insertions without this check use 'erase'
  void erase_existing(const T id_or_index) {
    assert(this->count(id_or_index) == 1);
    this->erase(id_or_index);
  }

  void dump(const std::string comment = "") const {
    std::cout << comment << ": ";
    int cnt = 0;
    for (const T& idx: *this) {
      std::cout << idx << ", ";

      if (cnt %20 == 0 && cnt != 0) {
        std::cout << "\n";
      }
      cnt++;
    }
    std::cout << "\n";
  }
};

#endif // __SHARED_DEFINES_H__
