// macros for logging must be in a header file.
// this has to be included after the logger module is imported.

#define ERROR(mess) print(1,mess,__FILE__,__func__,__LINE__)
#define WARN(mess) print(2,mess,__FILE__,__func__,__LINE__)
#define INFO(mess) print(3,mess,__FILE__,__func__,__LINE__)
#define DEBUG(mess) print(4,mess,__FILE__,__func__,__LINE__)
#define TRACE(mess) print(5,mess,__FILE__,__func__,__LINE__)
