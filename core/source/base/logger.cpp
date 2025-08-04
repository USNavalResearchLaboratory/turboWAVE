module;

#include "tw_includes.h"

/// module logger works together with macros in tw_logger.h.
/// usage example: `logger::DEBUG(std::format("result was {}",result))`
/// level is controlled by `TW_LOG` environment variable.
/// levels are error,warn,info,debug,trace.
export module logger;
import base;
int tw_log_level = 1; // 0=none,1=err,2=warn,3=info,4=debug,5=trace

export namespace logger {

    void init() {
        const char * raw = std::getenv("TW_LOG");
        std::string twenv;
        if (!raw) {
            twenv = "error";
        } else {
            twenv = std::string(raw);
        }
        if (twenv=="error") {
            tw_log_level = 1;
        } else if (twenv=="warn") {
            tw_log_level = 2;
        } else if (twenv=="info") {
            tw_log_level = 3;
        } else if (twenv=="debug") {
            tw_log_level = 4;
        } else if (twenv=="trace") {
            tw_log_level = 5;
        } else {
            tw_log_level = 0;
        }
    }
    void print(int level,const std::string& mess,const std::string& path,const std::string& func,int line) {
        if (tw_log_level < level) {
            return;
        }
        std::string file = path.substr(path.rfind('/')+1);
        switch (level) {
            case 1:
                std::println(std::cout,"{}[{}][{}][{}]: {}",
                    term::error,file,func,line,mess);
                break;
            case 2:
                std::println(std::cout,"{}[{}][{}][{}]: {}",
                    term::warning,file,func,line,mess);
                break;
            case 3:
                std::println(std::cout,"{}INFO{}[{}][{}][{}]: {}",
                    term::green,term::reset_all,file,func,line,mess);
                break;
            case 4:
                std::println(std::cout,"{}DEBUG{}[{}][{}][{}]: {}",
                    term::blue,term::reset_all,file,func,line,mess);
                break;
            case 5:
                std::println(std::cout,"{}TRACE{}[{}][{}][{}]: {}",
                    term::cyan,term::reset_all,file,func,line,mess);
                break;

        }
        std::flush(std::cout);
    }
}
