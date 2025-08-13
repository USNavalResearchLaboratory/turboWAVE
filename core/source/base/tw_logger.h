// Macros for logging.
// When expanded these call into the `base/logger` module.

#define ERROR(mess) print(1,mess,__FILE__,__func__,__LINE__)
#define WARN(mess) print(2,mess,__FILE__,__func__,__LINE__)
#define INFO(mess) print(3,mess,__FILE__,__func__,__LINE__)
#define DEBUG(mess) print(4,mess,__FILE__,__func__,__LINE__)
#define TRACE(mess) print(5,mess,__FILE__,__func__,__LINE__)

#ifdef ENABLE_BOUNDS_CHECKING
#define BOUNDS(idx)\
if (idx > array.size() || idx < 0) {\
    throw tw::FatalError(std::format("access violation {}/{}",idx,array.size()));\
}
#else
#define BOUNDS(idx)
#endif