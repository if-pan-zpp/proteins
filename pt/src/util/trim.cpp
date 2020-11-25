#include "trim.h"
#include <algorithm>
using namespace std;

namespace util {
    static bool notspace(unsigned char ch) {
        return !std::isspace(ch);
    }
}

std::string util::ltrim(std::string s) {
    if (s.empty()) return s;
    s.erase(s.begin(), find_if(s.begin(), s.end(), util::notspace));
    return s;
}

std::string util::rtrim(std::string s) {
    if (s.empty()) return s;
    s.erase(find_if(s.rbegin(), s.rend(), util::notspace).base(), s.end());
    return s;
}

std::string util::trim(std::string s) {
    return ltrim(rtrim(move(s)));
}
