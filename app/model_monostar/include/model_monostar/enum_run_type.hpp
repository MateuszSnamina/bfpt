#pragma once

// STD:
#include <string>
#include <ostream>
#include <stdexcept>

enum class RunType {
    G,
    E,
    EG
};

inline std::ostream& operator<<(std::ostream& s, const RunType& r) {
    if (r == RunType::G) {
        s << "RunType::G";
    } else if (r == RunType::E) {
        s << "ModelType::E";
    } else if (r == RunType::EG) {
        s << "ModelType::EG";
    } else {
        throw std::logic_error("Invalid run type enum value!");
    }
    return s;
}
