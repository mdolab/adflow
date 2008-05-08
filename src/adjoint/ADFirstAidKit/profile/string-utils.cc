#include "string-utils.h"

using namespace std;

string prefix (const string & str, int len)
{
        if (len >= 0)
                return str.substr (0, len);
        else
                return str.substr (0, str.length () + len);
}

string suffix (const string & str, int len)
{
        if (len >= 0)
                return str.substr (str.length () - len, len);
        else
                return str.substr (-len, str.length () + len);
}

vector < string > split (const string & str, char separator)
{
        unsigned int components = 0;
        unsigned int i = 0;
        unsigned int l = str.length ();
        // Count the components
        while (true) {
                while (i < l && str[i] == separator)
                        i++;
                if (i >= l)
                        break;
                components++;
                while (i < l && str[i] != separator)
                        i++;
        }
        // Allocate the vector
        vector < string > split_list;
	split_list.reserve(4);
        // Split!
        int lpos = 0, pos = 0;
        while (true) {
                while (pos < l && str[pos] == separator)
                        pos++;
                if (pos >= l)
                        break;
                lpos = pos;
                while (pos < l && str[pos] != separator)
                        pos++;
                split_list.push_back (str.substr (lpos, pos - lpos));
        }
        return split_list;
}
string last_part_after(const string &str, char separator)
{
        int pos = str.find_last_of(separator);
        return str.substr(pos+1);
}
