#pragma once

// STD:
#include <string>
#include <ostream>
#include <stdexcept>

namespace so_app {

enum class ModelType {
    AFFO
};

inline std::ostream& operator<<(std::ostream& s, const ModelType& m) {
    if (m == ModelType::AFFO) {
        s << "ModelType::AFFO";
    } else {
        throw std::logic_error("Invalid model enum value!");
    }
    return s;
}

}  // end of namespace so_app
