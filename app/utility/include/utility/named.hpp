#ifndef UTILITY_NAMED_HPP
#define UTILITY_NAMED_HPP

#include<string>

namespace  utility {

template<typename T>
struct Named {
    Named(T value, std::string name = "unnamed") :
    name(name),
    value(value) {
    }
    std::string name;
    T value;
};

}

#endif // UTILITY_NAMED_HPP
