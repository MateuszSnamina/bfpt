#pragma once

// STD:
#include <string>
#include <ostream>
#include <stdexcept>

enum class ModelType {
    FM,
    AF,
    FO
};

inline std::ostream& operator<<(std::ostream& s, const ModelType& m) {
    if (m == ModelType::FM) {
        s << "ModelType::FM";
    } else if (m == ModelType::AF) {
        s << "ModelType::AF";
    } else if (m == ModelType::AF) {
        s << "ModelType::FO";
    } else {
        throw std::logic_error("Invalid model enum value!");
    }
    return s;
}
