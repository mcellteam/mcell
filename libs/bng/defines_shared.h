/**
 * Definitions shared among BNG and other sources that might use it.
 *
 * No namespace is used here (except for float_t that in enclosed in namespace Common).
 */

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
#include <boost/container/flat_map.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

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

#include "bng_config.h"

const Common::float_t BNG_PI = 3.14159265358979323846;
const Common::float_t BNG_N_AV = 6.0221417930e23;


#if FLOAT_T_BYTES == 8
const Common::float_t EPS = 1e-12; // same as EPS_C
const Common::float_t SQRT_EPS = 1e-6;
const Common::float_t GIGANTIC4 = 1e140;
#else
#error "Base type float32 is not supported yet"
#endif

// ---------------------------------- fixed constants and specific typedefs -------------------
const Common::float_t POS_INVALID = FLT_MAX; // cannot be NAN because we cannot do any comparison with NANs

const Common::float_t FLT_INVALID = FLT_MAX;

const Common::float_t TIME_INVALID = -1;
const Common::float_t TIME_FOREVER = FLT_MAX; // this max is sufficient for both float and double
const Common::float_t TIME_SIMULATION_START = 0;

const Common::float_t UINT_INVALID = UINT32_MAX; // invalid value to be used for any invalid unsigned integer values
const Common::float_t UINT_INVALID2 = UINT32_MAX - 1; // second invalid value not to be used, in this case for any purpose

const uint ID_INVALID = UINT_INVALID; // general invalid index, should not be used when a definition for a specific type is available
const uint INDEX_INVALID = UINT_INVALID; // general invalid index, should not be used when a definition for a specific type is available


// molecule id is a unique identifier of a molecule,
// no 2 molecules may have the same ID in the course of a simulation
// not used from the BNG lib yet
typedef uint molecule_id_t;
const molecule_id_t MOLECULE_ID_INVALID = ID_INVALID;

// unique species id,
// every distinct species that exists or
typedef uint species_id_t;
const species_id_t SPECIES_ID_INVALID = ID_INVALID;

// unique pattern id
typedef uint pattern_id_t;
const pattern_id_t PATTERN_ID_INVALID = ID_INVALID;

typedef int32_t orientation_t;
const orientation_t ORIENTATION_DOWN = -1;
const orientation_t ORIENTATION_NONE = 0;
const orientation_t ORIENTATION_UP = 1;
const orientation_t ORIENTATION_NOT_SET = 2;

// index of reaction in a reaction class
typedef int reaction_index_t;


#ifndef NDEBUG
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
 * Set based on boost's flat_set. 
 * Used to hold sets of ids or indices of molecules or other items,
 * extended to check for unique insertions and checked removals.
 *
 * T is really not needed since it is an uint, but the declaration then
 * shows what type of values are expected.
 */
template<typename Key>
class uint_set: public base_flat_set<Key> {
public:
  // insert with check that the item is not there yet
  // for insertions without this check use 'insert'
  void insert_unique(const Key id_or_index) {
    assert(this->count(id_or_index) == 0);
    this->insert(id_or_index);
  }

  // erase with check that the item is present
  // for insertions without this check use 'erase'
  void erase_existing(const Key id_or_index) {
    assert(this->count(id_or_index) == 1);
    this->erase(id_or_index);
  }

  void dump(const std::string comment = "") const;
};


/**
 * Set based on googles's dense_hash_set.
 */
template<typename Key>
class uint_dense_hash_set: public google::dense_hash_set<Key> {
public:
  uint_dense_hash_set() {
    // dense_hash_set requires this to be called
    this->set_empty_key(UINT_INVALID);
    this->set_deleted_key(UINT_INVALID2);
  }

  // insert with check that the item is not there yet
  // for insertions without this check use 'insert'
  void insert_unique(const Key id_or_index) {
    assert(this->count(id_or_index) == 0);
    this->insert(id_or_index);
  }

  // erase with check that the item is present
  // for insertions without this check use 'erase'
  void erase_existing(const Key id_or_index) {
    assert(this->count(id_or_index) == 1);
    this->erase(id_or_index);
  }
};


/**
 * Map based on boost's flat_map.
 * 
 * Not used yet
 */ 
template<typename Key, typename Value>
class uint_flat_map: public boost::container::flat_map<Key, Value> {

};

/**
 * Map based on googles's dense_hash_map.
 * 
 *  Warning: value must be POD (plain old data), this is not checked during compilations.
 *
 *  Not used yet
 */
template<typename Key, typename Value>
class uint_dense_hash_map: public google::dense_hash_map<Key, Value> {
public:
  uint_dense_hash_map() {
    assert(sizeof(Value) <= 8);
    // dense_hash_map requires this to be called
    this->set_empty_key(UINT_INVALID);
    this->set_deleted_key(UINT_INVALID2);
  }
};



template<typename Key>
void uint_set<Key>::dump(const std::string comment) const {
  std::cout << comment << ": ";
  int cnt = 0;
  for (const Key& idx: *this) {
    std::cout << idx << ", ";

    if (cnt %20 == 0 && cnt != 0) {
      std::cout << "\n";
    }
    cnt++;
  }
  std::cout << "\n";
}


#endif // __SHARED_DEFINES_H__
