#pragma once

#include <kbasis/basis.hpp>

#include <kstate_impl/kstate_streamer.hpp>

#include <utility/is_base_of_template.hpp>
#include <utility/remove_cvref.hpp>

#include <boost/range/adaptor/indirected.hpp>

#include <iomanip>

// #######################################################################
// ##  BasisStreamer                                                    ##
// #######################################################################

namespace kbasis {

// ***********************************************************************

template<typename _BasisT>
class BasisStreamer;

// ***********************************************************************

template <typename BasisT>
BasisStreamer<BasisT>
make_basis_streamer(
        BasisT&&,
        const extension::boost::RangeStreamerSettings<typename utility::remove_cvref_t<BasisT>::KstateT>);

// ***********************************************************************

template<typename _BasisT>
class BasisStreamer {
    static_assert(!std::is_array_v<_BasisT>);
    static_assert(!std::is_function_v<_BasisT>);
    static_assert(!std::is_void_v<std::decay<_BasisT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_BasisT>>);
    static_assert(!std::is_enum_v<std::decay<_BasisT>>);
    static_assert(!std::is_union_v<std::decay<_BasisT>>);
    static_assert(std::is_class_v<std::decay<_BasisT>>);
    static_assert(!std::is_pointer_v<std::decay<_BasisT>>);
    static_assert(!std::is_member_object_pointer_v<_BasisT>);
    static_assert(!std::is_member_function_pointer_v<_BasisT>);
    static_assert(!std::is_volatile_v<_BasisT>);
    static_assert(!std::is_rvalue_reference_v<_BasisT>,
    "BasisT must not be a rvalue reference.");
    static_assert(std::is_lvalue_reference_v<_BasisT> || (!std::is_reference_v<_BasisT> && !std::is_const_v<_BasisT>),
    "BasisT must be of form: `T`, `T&` or `const T&`.");
    static_assert(utility::is_base_of_template_v<utility::remove_cvref_t<_BasisT>, Basis>,
    "remove_cvref_t<BasisT> must be derived from Basis");
public: // Helper types:
    using BasisT = _BasisT;
    using KstateT = typename utility::remove_cvref_t<BasisT>::KstateT;
public: // Factory function:
    friend
    BasisStreamer<BasisT> make_basis_streamer<BasisT>(BasisT&&, const extension::boost::RangeStreamerSettings<KstateT>);
public: // API:
    ::std::ostream& stream(std::ostream&) const;
    std::string str() const;
private:
    BasisStreamer(std::add_lvalue_reference_t<BasisT>, extension::boost::RangeStreamerSettings<KstateT>);
private:
    const _BasisT _basis;
    const extension::boost::RangeStreamerSettings<KstateT> _range_streamer_settings;
};

// ***********************************************************************

template<typename _BasisT>
BasisStreamer<_BasisT>::BasisStreamer(
        std::add_lvalue_reference_t<_BasisT> basis,
        extension::boost::RangeStreamerSettings<KstateT> range_streamer_settings) :
    _basis(std::forward<_BasisT>(basis)),
    _range_streamer_settings(range_streamer_settings) {
}

// ***********************************************************************

template<typename _BasisT>
::std::ostream&
BasisStreamer<_BasisT>::stream(std::ostream& os) const {
    // Defaults for basis range streamer:
    const ::std::function<void(::std::ostream&)> default_stream_preparer =
            [](std::ostream& s) { s << "𝔹𝔸𝕊𝕀𝕊-BEGIN" << std::endl; };
    const ::std::function<void(::std::ostream&, size_t)> default_stream_sustainer =
            [](std::ostream& s, size_t i) { s << " - " << std::right << std::setw(8) << i << " : "; };
    const ::std::function<void(::std::ostream&, KstateT)> default_stream_value_putter =
            [](::std::ostream& s, KstateT) { s << "[PLACE-FOR-KSTATE]"; };
    const ::std::function<void(::std::ostream&)> default_stream_separer =
            [](std::ostream& s) { s << std::endl; };
    const ::std::function<void(::std::ostream&)> default_stream_finisher =
            [](std::ostream& s) { s << std::endl << "𝔹𝔸𝕊𝕀𝕊-END" << std::endl; };
    bool default_format_independence_flag = true;
    // Apply overrules:
    const auto stream_preparer = (
                _range_streamer_settings._stream_preparer ?
                    *_range_streamer_settings._stream_preparer :
                    default_stream_preparer);
    const auto stream_sustainer = (
                _range_streamer_settings._stream_sustainer ?
                    *_range_streamer_settings._stream_sustainer :
                    default_stream_sustainer);
    const auto stream_value_putter = (
                _range_streamer_settings._stream_value_putter ?
                    *_range_streamer_settings._stream_value_putter :
                    default_stream_value_putter);
    const auto stream_separer = (
                _range_streamer_settings._stream_separer ?
                    *_range_streamer_settings._stream_separer :
                    default_stream_separer);
    const auto stream_finisher = (
                _range_streamer_settings._stream_finisher ?
                    *_range_streamer_settings._stream_finisher :
                    default_stream_finisher);
    const auto format_independence_flag = (
                _range_streamer_settings._format_independence_flag ?
                    *_range_streamer_settings._format_independence_flag :
                    default_format_independence_flag);
    // Stream:
    extension::boost::stream_range_impl<decltype(_basis.vec_index() | ::boost::adaptors::indirected)>(
                _basis.vec_index() | ::boost::adaptors::indirected,
                os,
                stream_preparer,
                stream_sustainer,
                stream_value_putter,
                stream_separer,
                stream_finisher,
                format_independence_flag);
    // Return:
    return os;
}

// ***********************************************************************

template<typename _BasisT>
std::string
BasisStreamer<_BasisT>::str() const {
    std::ostringstream oss;
    stream(oss);
    return oss.str();
}

// ***********************************************************************

template<typename BasisT>
BasisStreamer<BasisT> make_basis_streamer(
        BasisT&& basis,
        const extension::boost::RangeStreamerSettings<typename utility::remove_cvref_t<BasisT>::KstateT> range_streamer_settings) {
    return BasisStreamer<BasisT>(std::forward<BasisT>(basis), range_streamer_settings);
}

} // end of namespace kbasis

// #######################################################################
// ##  pragma: operator||, operator<<                                   ##
// #######################################################################

namespace kbasis::pramga {

template<typename BasisT>
BasisStreamer<BasisT>
operator&&(
        BasisT&& basis,
        const extension::boost::RangeStreamerSettings<typename utility::remove_cvref_t<BasisT>::KstateT> range_streamer_settings) {
    return make_basis_streamer<BasisT>(std::forward<BasisT>(basis), range_streamer_settings);
}

template<typename BasisT>
std::ostream& operator<<(
        std::ostream& os,
        const BasisStreamer<BasisT>& basis_streamer) {
    basis_streamer.stream(os);
    return os;
}

}  // namespace kbasis::pramga
