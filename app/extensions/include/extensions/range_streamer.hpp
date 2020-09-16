#pragma once

#include <extensions/stream_fromat_stacker.hpp>

#include <boost/range.hpp>
#include <boost/range/adaptor/indexed.hpp>

#include <functional>
#include <optional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <type_traits>

// #######################################################################
// ##  RangeStreamerSettings                                            ##
// #######################################################################

namespace extension::boost {

template<typename _T>
struct RangeStreamerSettings {
    // Helper types:
    using T = _T;
    // Setters - general:
    RangeStreamerSettings<T>& set_stream_preparer(::std::function<void(::std::ostream&)>);
    RangeStreamerSettings<T>& set_stream_sustainer(::std::function<void(::std::ostream&, size_t)>);
    RangeStreamerSettings<T>& set_stream_value_putter(::std::function<void(::std::ostream&, T)>);
    RangeStreamerSettings<T>& set_stream_separer(::std::function<void(::std::ostream&)>);
    RangeStreamerSettings<T>& set_stream_finisher(::std::function<void(::std::ostream&)>);
    RangeStreamerSettings<T>& set_format_independence_flag(bool = true);
    // Setters - convenience:
    RangeStreamerSettings<T>& set_null_preparer();
    RangeStreamerSettings<T>& set_char_preparer(char c);
    RangeStreamerSettings<T>& set_string_preparer(::std::string str);
    RangeStreamerSettings<T>& set_null_sustainer();
    RangeStreamerSettings<T>& set_null_separer();
    RangeStreamerSettings<T>& set_char_separer(char c);
    RangeStreamerSettings<T>& set_string_separer(::std::string str);
    RangeStreamerSettings<T>& set_null_finisher();
    RangeStreamerSettings<T>& set_char_finisher(char c);
    RangeStreamerSettings<T>& set_string_finisher(::std::string str);
    RangeStreamerSettings<T>& in_bracket_round();
    RangeStreamerSettings<T>& in_bracket_curly();
    RangeStreamerSettings<T>& in_bracket_square();
    RangeStreamerSettings<T>& in_bracket_angle();
    RangeStreamerSettings<T>& in_quotation_single();
    RangeStreamerSettings<T>& in_quotation_double();
    RangeStreamerSettings<T>& in_quotation_back();
    RangeStreamerSettings<T>& in_null();
    RangeStreamerSettings<T>& like_python_tuple();
    RangeStreamerSettings<T>& like_python_list();
    RangeStreamerSettings<T>& like_python_set();
    // Fields:
    ::std::optional<::std::function<void(::std::ostream&)>> _stream_preparer;
    ::std::optional<::std::function<void(::std::ostream&, size_t)>> _stream_sustainer;
    ::std::optional<::std::function<void(::std::ostream&, T)>> _stream_value_putter;
    ::std::optional<::std::function<void(::std::ostream&)>> _stream_separer;
    ::std::optional<::std::function<void(::std::ostream&)>> _stream_finisher;
    ::std::optional<bool> _format_independence_flag;
};

// ***********************************************************************

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_preparer(::std::function<void(::std::ostream&)> _) {
    _stream_preparer = _;
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_sustainer(::std::function<void(::std::ostream&, size_t)> _) {
    _stream_sustainer = _;
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_value_putter(::std::function<void(::std::ostream&, T)> _) {
    _stream_value_putter = _;
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_separer(::std::function<void(::std::ostream&)> _) {
    _stream_separer = _;
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_format_independence_flag(bool _) {
    _format_independence_flag = _;
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_finisher(::std::function<void(::std::ostream&)> _) {
    _stream_finisher = _;
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_null_preparer() {
    _stream_preparer = [](::std::ostream&) {};
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_char_preparer(char c) {
    _stream_preparer = [c](::std::ostream& s) { s << c; };
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_string_preparer(::std::string str) {
    _stream_preparer = [str](::std::ostream& s) { s << str; };
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_null_sustainer() {
    _stream_sustainer = [](::std::ostream&, size_t) {};
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_null_separer() {
    _stream_separer = [](::std::ostream&) {};
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_char_separer(char c) {
    _stream_separer = [c](::std::ostream& s) { s << c; };
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_string_separer(::std::string str) {
    _stream_separer = [str](::std::ostream& s) { s << str; };
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_null_finisher() {
    _stream_finisher = [](::std::ostream&) {};
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_char_finisher(char c) {
    _stream_finisher = [c](::std::ostream& s) { s << c; };
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_string_finisher(::std::string str) {
    _stream_finisher = [str](::std::ostream& s) { s << str; };
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_bracket_round() {
    set_char_preparer('(');
    set_char_finisher(')');
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_bracket_curly() {
    set_char_preparer('{');
    set_char_finisher('}');
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_bracket_square() {
    set_char_preparer('[');
    set_char_finisher(']');
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_bracket_angle() {
    set_char_preparer('<');
    set_char_finisher('>');
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_quotation_single() {
    set_char_preparer('\'');
    set_char_finisher('\'');
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_quotation_double() {
    set_char_preparer('"');
    set_char_finisher('"');
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_quotation_back() {
    set_char_preparer('`');
    set_char_finisher('`');
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_null() {
    set_char_preparer(' ');
    set_char_finisher(' ');
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::like_python_tuple() {
    set_char_preparer('[');
    set_char_finisher(']');
    set_string_separer(", ");
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::like_python_list() {
    set_char_preparer('[');
    set_char_finisher(']');
    set_string_separer(", ");
    return *this;
}

template<typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::like_python_set() {
    set_char_preparer('{');
    set_char_finisher('}');
    set_string_separer(", ");
    return *this;
}

}  // namespace extension::boost

// #######################################################################
// ##  stream_range_impl                                                ##
// #######################################################################

namespace extension::boost {

template<typename R>
void stream_range_impl(
        const ::std::remove_reference_t<R>& range,
        ::std::ostream& os,
        ::std::function<void(::std::ostream&)> stream_preparer,
        ::std::function<void(::std::ostream&, size_t)> stream_sustainer,
        ::std::function<void(::std::ostream&, typename ::boost::range_value<R>::type)> stream_value_putter,
        ::std::function<void(::std::ostream&)> stream_separer,
        ::std::function<void(::std::ostream&)> stream_finisher,
        bool format_independence_flag) {
    const extension::std::StreamFromatStacker stream_format_stacker(os, format_independence_flag);
    stream_preparer(os);
    for (const auto& _ : range | ::boost::adaptors::indexed(0)) {
        if (_.index() != 0) {
            stream_separer(os);
        }
        stream_sustainer(os, _.index());
        stream_value_putter(os, _.value());
    }
    stream_finisher(os);
}

} // namespace extension::boost

// #######################################################################
// ##  RangeStreamer                                                    ##
// #######################################################################

namespace extension::boost {

template<class _R>
class RangeStreamer;

template<class R>
RangeStreamer<R> make_range_streamer(
        R&& range,
        const RangeStreamerSettings<typename ::boost::range_value<R>::type> range_streamer_settings);

template<class _R>
class RangeStreamer {
    static_assert(! ::std::is_rvalue_reference_v<_R>,
    "R must not be a rvalue reference.");
    static_assert(::std::is_lvalue_reference_v<_R> || (!::std::is_reference_v<_R> && !::std::is_const_v<_R>),
    "R must be of form: `T`, `T&` or `const T&`.");
public: // Helper types:
    using R = _R;
    using T = typename ::boost::range_value<R>::type;
public: // Factory function:
    friend RangeStreamer<R> make_range_streamer<R>(R&& range, const RangeStreamerSettings<T> range_streamer_settings);
public: // API:
    ::std::ostream& stream(::std::ostream& os) const;
    ::std::string str() const;
private:
    RangeStreamer(::std::add_rvalue_reference_t<R> range, const RangeStreamerSettings<T> range_streamer_settings);
private:
    const R _range;
    const RangeStreamerSettings<typename ::boost::range_value<R>::type> _range_streamer_settings;
};

template<class _R>
RangeStreamer<_R>::RangeStreamer(
        ::std::add_rvalue_reference_t<R> range,
        const RangeStreamerSettings<T> range_streamer_settings) :
    _range(::std::forward<R>(range)),
    _range_streamer_settings(range_streamer_settings) {
}


template<class _R>
::std::ostream&
RangeStreamer<_R>::stream(::std::ostream& os) const {
    // Defaults:
    const auto default_stream_preparer = [](::std::ostream& s) { s << "{"; };
    const auto default_stream_sustainer = [](::std::ostream& s, size_t i) { s << i << ":"; };
    const auto default_stream_value_putter = [](::std::ostream& s, T t) { s << t; };
    const auto default_stream_separer = [](::std::ostream& s) { s << ","; };
    const auto default_stream_finisher = [](::std::ostream& s) { s << "}"; };
    const bool default_format_independence_flag = true;
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
    stream_range_impl<R>(
                _range,
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

template<class _R>
::std::string RangeStreamer<_R>::str() const {
    ::std::ostringstream oss;
    stream(oss);
    return oss.str();
}


template<class R>
RangeStreamer<R> make_range_streamer(R&& range, const RangeStreamerSettings<typename ::boost::range_value<R>::type> range_streamer_settings) {
    return RangeStreamer<R>(::std::forward<R>(range), range_streamer_settings);
}

}  // namespace extension::boost

// #######################################################################
// ##  operator|, Stringifier                                           ##
// #######################################################################

namespace extension::boost {
namespace stream_pragma {

template<class R>
RangeStreamer<R> operator|(R&& range, const RangeStreamerSettings<typename ::boost::range_value<R>::type> range_streamer_settings){
    return make_range_streamer(::std::forward<R>(range), range_streamer_settings);
}

template<class R>
::std::ostream& operator<<(::std::ostream& os, const RangeStreamer<R>& range_streamer) {
    return range_streamer.stream(os);
}

struct Stringifier {
};

template<class R>
::std::string operator*(const RangeStreamer<R>& range_streamer, Stringifier) {
    return range_streamer.str();
}

template<typename T>
using RSS = RangeStreamerSettings<T>;

}  // namespace stream_pragma
}  // namespace extension::boost
