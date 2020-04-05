#ifndef SRC4_COUNTED_VOLUME_UTIL_H_
#define SRC4_COUNTED_VOLUME_UTIL_H_

namespace MCell {

class World;

namespace CountedVolumesUtil {

// return true if counted volumes were correctly set up
bool initialize_counted_volumes(World* world);

};

} // namespace MCell

#endif // SRC4_COUNTED_VOLUME_UTIL_H_
