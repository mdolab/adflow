#include "reader.h"
#include "string-utils.h"
#include <cstdio>
#include <stdexcept>
#include <map>

// Real impl
using namespace std;
struct TapenadeLogReaderImpl {
        FILE *file;
	int file_length;
	map<int,string> function_id_2_name;
	~TapenadeLogReaderImpl() {
		if (file)
			fclose(file);
	}
};
template <typename T> static T read(TapenadeLogReaderImpl *impl) {
  T data;
  fread(&data, 1, sizeof(T), impl->file);
  return data;
}
template <> static inline string read(TapenadeLogReaderImpl *impl) {
  unsigned int len = read<unsigned int>(impl);
  char *buffer = new char[len];
  fread(buffer, 1, len, impl->file);
  string ret(buffer, len);
  delete[] buffer;
  return ret;
}
template <> static inline TapenadeLogEvent read(TapenadeLogReaderImpl *impl) {
	unsigned int function_id = read<unsigned int>(impl);
	string &function_name = impl->function_id_2_name[function_id];
	unsigned int kind, size;
	unsigned long long int time;
	kind = read<unsigned int>(impl);
	time = read<unsigned long long int>(impl);
	size = read<unsigned int>(impl);
	return TapenadeLogEvent(function_name, kind, time, size);
}
static inline bool eof(TapenadeLogReaderImpl *impl) {
	int curPos = ftell(impl->file);
	return (curPos >= impl->file_length);
}

static void initialize(TapenadeLogReaderImpl *&impl, const string &filename, bool &ok, bool &init) throw(runtime_error) {
        if ((! init) && ok) {
                impl = new TapenadeLogReaderImpl();
                impl->file = fopen(filename.c_str(), "r");
                if (! impl->file) {
                  ok = false;
                  throw new runtime_error("Can't open "+
		                           filename+
		                          "for reading");
		}
		fseek(impl->file, 0, SEEK_END);
		impl->file_length = ftell(impl->file);
		rewind(impl->file);
		// Read functions ids
		unsigned int hashtable_count;
		hashtable_count = read<unsigned int>(impl);
		while (hashtable_count--) {
			unsigned int function_id = read<unsigned int>(impl);
			impl->function_id_2_name[function_id] =
			  read<string>(impl);
		}
		init = true;
	}
}

// Interface
#define TapenadeLogReaderTemplate(type, name) \
type TapenadeLogReader::name() throw(runtime_error) {\
        if (! initialized)\
                initialize(impl, filename, ok, initialized);\
        if (ok)\
                return read<type>(impl);\
	else \
	        throw runtime_error("End of file reached in \""+\
		                    filename+"\"");\
}
bool TapenadeLogReader::eof() throw(runtime_error) {
        if (! initialized)
                initialize(impl, filename, ok, initialized);
        if (ok)
                return ::eof(impl);
	else 
	        throw runtime_error("End of file reached in \""+
		                    filename+"\"");
}
TapenadeLogReader::~TapenadeLogReader() {
	if (impl) {
		delete impl;
	}
}

TapenadeLogReaderTemplate(TapenadeLogEvent, readEvent);

// Interface for ReadingOnlyInAdjointBranch
ReadingOnlyInAdjointBranch::ReadingOnlyInAdjointBranch(const std::string &filename)
{
        unread = false;
        done =false;
        adjoint_branch = false;
        reader = new TapenadeLogReader(filename);
}
ReadingOnlyInAdjointBranch::~ReadingOnlyInAdjointBranch() {
        if (unreadEvent)
                delete unreadEvent;
        if (reader)
                delete reader;
}

bool ReadingOnlyInAdjointBranch::eof() throw(runtime_error)
{
        if (adjoint_branch)
                return reader->eof();
        if (done)
                return true;

        while (! reader->eof()) {
                TapenadeLogEvent event(reader->readEvent());
                string suffix = last_part_after(event.function_name, '_');
                if (suffix == "b" || suffix == "fwd") {
                        adjoint_branch_first_function_name = event.function_name;
                        adjoint_branch_first_function_recursion = 1;
                        unreadEvent = new TapenadeLogEvent(event);
                        unread = true;
                        adjoint_branch = true;
                        return false;
                }
        }
        return true;
}
TapenadeLogEvent ReadingOnlyInAdjointBranch::readEvent() throw(runtime_error)
{
        if (! adjoint_branch)
                throw runtime_error("readEvent called without verifying eof!");
        if (unread) {
                TapenadeLogEvent event(*unreadEvent);
                delete unreadEvent;
                unreadEvent = 0;
                unread = false;
                return event;
        }
        TapenadeLogEvent event(reader->readEvent());
        if (event.function_name == adjoint_branch_first_function_name) {
                if (event.kind == END)
                        if (adjoint_branch_first_function_recursion == 1) {
                                adjoint_branch = false;
                                done = true;
                        } else
                                adjoint_branch_first_function_recursion--;
                if (event.kind == BEGIN)
                        adjoint_branch_first_function_recursion++;
        }
        return event;
}
