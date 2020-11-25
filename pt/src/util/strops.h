#pragma once
#include <string>
#include <vector>

namespace util {
    std::string_view ltrim(std::string_view s);
    std::string_view rtrim(std::string_view s);
    std::string_view trim(std::string_view s);
    std::vector<std::string_view> split(std::string_view s);
}
