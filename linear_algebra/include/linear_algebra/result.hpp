#ifndef LINEAR_ALGEBRA_RESULT_HPP
#define LINEAR_ALGEBRA_RESULT_HPP

#include <any>
#include <exception>
#include <functional>
#include <type_traits>
#include <variant>

#include <cassert>

// #######################################################################
// ## LinearAlgebraRuntimeExceptions                                    ##
// #######################################################################

namespace lin_alg {

/*
 * The type used in conjunction with Result.
 */
class LinearAlgebraRuntimeException : public std::runtime_error {
   public:
    LinearAlgebraRuntimeException(const std::string& message, std::any datails);
    std::any get_details() const;

   private:
    const std::any _datails;
};

inline LinearAlgebraRuntimeException::LinearAlgebraRuntimeException(const std::string& message, std::any datails)
    : runtime_error(message),
      _datails(datails) {
}

inline std::any LinearAlgebraRuntimeException::get_details() const {
    return _datails;
}

}  // namespace lin_alg

// #######################################################################
// ## Result                                                            ##
// #######################################################################

namespace lin_alg {

template <typename OutT, typename ExceptionT>
class Result {
    static_assert(std::is_base_of<std::exception, ExceptionT>::value);

   public:
    static Result<OutT, ExceptionT> Ok(OutT);
    static Result<OutT, ExceptionT> Err(ExceptionT);
    bool is_ok() const;
    bool is_err() const;
    OutT unwrap() const;
    OutT unwrap_or(OutT optb) const;
    OutT unwrap_or_else(std::function<OutT(ExceptionT)> op) const;

   private:
    Result(std::variant<OutT, ExceptionT>);
    const std::variant<OutT, ExceptionT> _variant;
};

}  // namespace lin_alg

// #######################################################################
// ## Result -- implementation                                          ##
// #######################################################################

namespace lin_alg {

// ------------------------------------------------------

template <typename OutT, typename ExceptionT>
Result<OutT, ExceptionT>::Result(std::variant<OutT, ExceptionT> variant)
    : _variant(variant) {
}

// ------------------------------------------------------

template <typename OutT, typename ExceptionT>
Result<OutT, ExceptionT> Result<OutT, ExceptionT>::Ok(OutT val) {
    const std::variant<OutT, ExceptionT> variant(std::in_place_type<OutT>, val);
    return Result(variant);
}

template <typename OutT, typename ExceptionT>
Result<OutT, ExceptionT> Result<OutT, ExceptionT>::Err(ExceptionT val) {
    const std::variant<OutT, ExceptionT> variant(std::in_place_type<ExceptionT>, val);
    return Result(variant);
}

// ------------------------------------------------------

template <typename OutT, typename ExceptionT>
bool Result<OutT, ExceptionT>::is_ok() const {
    return std::holds_alternative<OutT>(_variant);
}

template <typename OutT, typename ExceptionT>
bool Result<OutT, ExceptionT>::is_err() const {
    return std::holds_alternative<ExceptionT>(_variant);
}

// ------------------------------------------------------

template <typename OutT, typename ExceptionT>
struct UnwrapVisitor {
    OutT operator()(const OutT& value) const {
        return value;
    }
    OutT operator()(const ExceptionT& exception) const {
        throw exception;
    }
};

template <typename OutT, typename ExceptionT>
OutT Result<OutT, ExceptionT>::unwrap() const {
    UnwrapVisitor<OutT, ExceptionT> visitor;
    return std::visit(visitor, _variant);
}

// ------------------------------------------------------
template <typename OutT, typename ExceptionT>
class UnwrapOrVisitor {
   public:
    UnwrapOrVisitor(const OutT& optb) : _optb(optb) {
    }
    OutT operator()(const OutT& value) const {
        return value;
    }
    OutT operator()(const ExceptionT&) const {
        return _optb;
    }

   private:
    const OutT& _optb;
};

template <typename OutT, typename ExceptionT>
OutT Result<OutT, ExceptionT>::unwrap_or(OutT optb) const {
    UnwrapOrVisitor<OutT, ExceptionT> visitor(optb);
    return std::visit(visitor, _variant);
}

// ------------------------------------------------------

template <typename OutT, typename ExceptionT>
class UnwrapOrElseVisitor {
   public:
    UnwrapOrElseVisitor(std::function<OutT(ExceptionT)> op) : _op(op) {
    }
    OutT operator()(const OutT& value) const {
        return value;
    }
    OutT operator()(const ExceptionT& exception) const {
        return op(exception);
    }

   private:
    std::function<OutT(ExceptionT)> _op;
};

template <typename OutT, typename ExceptionT>
OutT Result<OutT, ExceptionT>::unwrap_or_else(std::function<OutT(ExceptionT)> op) const {
    UnwrapOrElseVisitor<OutT, ExceptionT> visitor(op);
    return std::visit(visitor, _variant);
}

// ------------------------------------------------------

}  // namespace lin_alg

#endif