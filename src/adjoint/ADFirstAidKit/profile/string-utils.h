#ifndef STRING_UTILS_H
#define STRING_UTILS_H 1
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

std::string prefix(const std::string &, int);
std::string suffix(const std::string &, int);
std::vector<std::string> split(const std::string &, char);
// Specialization of split(str, sep).back()
std::string last_part_after(const std::string &, char);

template <typename T> std::string toString(const T &ref)
{
        using namespace std;
        ostringstream string_buffer;
        string_buffer << ref;
        return string_buffer.str();
}
#endif
