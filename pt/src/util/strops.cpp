#include "strops.h"
#include <algorithm>
using namespace std;

namespace util {
    static bool notspace(unsigned char ch) {
        return !std::isspace(ch);
    }
}

std::string_view util::ltrim(std::string_view s) {
    auto beg = find_if(s.begin(), s.end(), util::notspace);
    auto end = s.end();
    return string_view(beg, distance(beg, end));
}

std::string_view util::rtrim(std::string_view s) {
    auto beg = s.begin();
    auto end = find_if(s.rbegin(), s.rend(), util::notspace).base();
    return string_view(beg, distance(beg, end));
}

std::string_view util::trim(std::string_view s) {
    return ltrim(rtrim(s));
}

std::vector<std::string_view> util::split(string_view s) {
    vector<string_view> res;
    int beg = -1;
    for (int i = 0; i <= s.size(); ++i) {
        if (beg < 0 && i < s.size() && !isspace(s[i])) {
            beg = i;
        }
        else if (i == s.size() || isspace(s[i])) {
            if (beg < i) res.emplace_back(&s[beg], i-beg);
            beg = -1;
        }
    }
    return res;
}
