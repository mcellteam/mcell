// TODO: better file name and guard macro name
#ifndef API_COMMON_H
#define API_COMMON_H

#include <ostream>
#include <sstream>
#include <exception>

#include "defines.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;
#include <string>
#include <vector>

#include "../generated/gen_constants.h"

namespace MCell {
namespace API {

// semantic check error codes
enum class SemRes {
  UNCHECKED, // not checked yet
  OK,        // nothing to do
  MESSAGE,   // only
  WARNING,   // print a warning
  ERROR,     // terminate on initialization of the whole model/or when added to model
             // (terminating right now would mean that the user can fix only one error at a time)
  FATAL_ERROR // terminate right now, cannot continue
};

//
const float_t FLT_UNSET = FLT_MAX;
const int INT_UNSET = INT32_MAX;
const Vec3 VEC3_UNSET(POS_INVALID);
const void* const PTR_UNSET = nullptr;
const std::string STR_UNSET = "unset";

// specialization for each type
static inline bool is_set(const float_t a) {
  return a != FLT_UNSET;
}
static inline bool is_set(const int a) {
  return a != INT32_MAX; // Python does not use unsigned integers
}
static inline bool is_set(const Vec3& a) {
  return a != VEC3_UNSET;
}
static inline bool is_set(const void* a) {
  return a != PTR_UNSET;
}
static inline bool is_set(const std::string& a) {
  return a != STR_UNSET && a != "";
}
template<typename T>
static inline bool is_set(const std::shared_ptr<T>& a) {
  return a.use_count() == 0;
}
template<typename T>
static inline bool is_set(const std::vector<T>& a) {
  return !a.empty();
}

// NOTE: for some reason a definition in defines.h is not found, although it works fine for Vec3
static inline std::ostream & operator<<(std::ostream &out, const IVec3& a) {
  out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
  return out;
}


template<typename T>
static inline std::string vec_ptr_to_str(const std::vector<T>& arr, const std::string ind="") {
  std::stringstream ss;
  ss << "[\n";
  for (size_t i = 0; i < arr.size(); i++) {
    ss << ind << " " << i << ":(" << arr[i]->to_str(ind) << ")";
    if (i + 1 != arr.size()) {
      ss << ", ";
    }
    ss << "\n";
  }
  ss << ind << "]";
  return ss.str();
}


template<typename T>
static inline std::string vec_nonptr_to_str(const std::vector<T>& arr, const std::string ind="") {
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < arr.size(); i++) {
    ss << i << ":" << arr[i];
    if (i + 1 != arr.size()) {
      ss << ", ";
    }
  }
  ss << "]";
  return ss.str();
}


typedef std::invalid_argument SemanticException;


// base class for all classes that hold the model input data
class BaseDataClass {
public:
  BaseDataClass()
    : class_name(STR_UNSET), name(STR_UNSET) {
  }
  virtual ~BaseDataClass(){
  }

  // we are storing class name for reporting
  std::string class_name;
  virtual void set_class_name(const std::string& class_name_) {
    class_name = class_name_;
  }
  virtual const std::string& get_class_name() const {
    return class_name;
  }

  // every object defined by the MCell API might have its name
  std::string name;
  virtual void set_name(const std::string& name_) {
    name = name_;
  }
  virtual const std::string& get_name() const {
    return name;
  }

  // this method is used to identify this particular object in error messages
  virtual std::string get_object_name() const {
    return get_class_name() + " '" + get_name() + "'";
  }

  // invalid implementation, to be overridden in actual derived classes
  virtual void ctor_postprocess() { };

  // invalid implementation, to be overridden in actual derived classes
  virtual SemRes check_semantics(std::ostream& out) const {
    out << "Semantic check for a derived class is not implemented.\n";
    return SemRes::FATAL_ERROR; /* not implemented*/
  };

  // calls virtual method, usually no need to override
  virtual SemRes check_semantics_cerr() const {
    return check_semantics(std::cerr);
  };

  // empty implementation, to be overridden in actual derived classes
  virtual std::string to_str(const std::string ind="") const {
    return "String dump for a derived class is not implemented.";
  }

  // calls virtual method, usually no need to override
  virtual void dump() const {
    std::cout << to_str() << "\n";
  }
};

} // namespace API
} // namespace MCell

#endif // API_COMMON_H
