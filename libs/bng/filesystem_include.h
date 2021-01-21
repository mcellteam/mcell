
#if __GNUC__ < 8 && !defined(_MSC_VER)
// gcc 6 & 7 have filesystem still under the experimental features
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif
