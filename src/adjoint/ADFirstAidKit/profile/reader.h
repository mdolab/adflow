#ifndef READER_H
#define READER_H 1
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <iterator>
#include "numeric-utils.h"

struct TapenadeLogReaderImpl;
enum TapenadeLogEventKind {
        BEGIN = 1,
        END = 2,
        ENDFWD = 3,
        BEGINSNAP = 4,
        ENDSNAP = 5,
        ENDORIG = 6
};

class TapenadeLogEvent
{
public:
        std::string &function_name;
        unsigned int kind;
        unsigned long long int time;
        unsigned int size;
        TapenadeLogEvent(std::string &_function_name ,unsigned int _kind,
                         unsigned long long int _time, unsigned int _size) :
                        function_name(_function_name),
                        kind(_kind),
                        time(_time),
                        size(_size)
        { }
        bool operator==(const TapenadeLogEvent &that)
        {
                return (function_name == that.function_name) &&
                       (kind = that.kind) &&
                       (time == that.time) &&
                       (size == that.size);
        }
        TapenadeLogEvent &operator=(const TapenadeLogEvent &other)
        {
                function_name = other.function_name;
                kind = other.kind;
                time = other.time;
                size = other.size;
        }
        std::string &get_function_name()
        {
                return function_name;
        }
        unsigned int get_kind()
        {
                return kind;
        }
        unsigned long long int get_time()
        {
                return time;
        }
        unsigned int get_size()
        {
                return size;
        }
};

inline std::ostream &operator<<(std::ostream &o, TapenadeLogEvent &event)
{
        using namespace std;
        char s = ' ';
        o << setw(14) << event.time << s;
        o << setw(9) << event.size << s;
        o << setw(2) << event.kind << s;
        o << event.function_name;
}

class TapenadeLogReader
{
        const std::string filename;
        bool ok;
        bool initialized;
        TapenadeLogReaderImpl *impl;
public:
        TapenadeLogReader(const std::string &_filename) : filename(_filename), ok(true), initialized(false)
        { }
        TapenadeLogEvent readEvent() throw (std::runtime_error);
        bool eof() throw (std::runtime_error);
        ~TapenadeLogReader();
};
class ReadingOnlyInAdjointBranch
{
        bool unread;
        bool done ;
        TapenadeLogEvent *unreadEvent;
        bool adjoint_branch ;
        std::string adjoint_branch_first_function_name;
        int adjoint_branch_first_function_recursion ;
        TapenadeLogReader *reader;
public:
        ReadingOnlyInAdjointBranch(const std::string &);
        ~ReadingOnlyInAdjointBranch();
        bool eof() throw(std::runtime_error);
        TapenadeLogEvent readEvent() throw(std::runtime_error);
};


#endif
