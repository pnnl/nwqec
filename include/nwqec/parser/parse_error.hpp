#pragma once

#include <stdexcept>
#include <string>

namespace NWQEC
{
    class ParseError : public std::runtime_error
    {
    public:
        ParseError(const std::string &message, int line, int column)
            : std::runtime_error("Parse error at line " + std::to_string(line) +
                                 ", column " + std::to_string(column) + ": " + message) {}
    };
}
