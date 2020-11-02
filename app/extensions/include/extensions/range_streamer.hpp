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

template <typename _T>
struct RangeStreamerSettings {
    // Helper types:
    using T = _T;
    using StreamPreparerFunT = ::std::function<void(::std::ostream&)>;
    using StreamSustainerFunT = ::std::function<void(::std::ostream&, size_t)>;
    using StreamValuePutterFunT = ::std::function<void(::std::ostream&, T)>;
    using StreamSeparerFunT = ::std::function<void(::std::ostream&)>;
    using StreamFinisherFunT = ::std::function<void(::std::ostream&)>;
    // Setters - general:
    RangeStreamerSettings<T>& set_stream_preparer(StreamPreparerFunT);
    RangeStreamerSettings<T>& set_stream_sustainer(StreamSustainerFunT);
    RangeStreamerSettings<T>& set_stream_value_putter(StreamValuePutterFunT);
    RangeStreamerSettings<T>& set_stream_separer(StreamSeparerFunT);
    RangeStreamerSettings<T>& set_stream_finisher(StreamFinisherFunT);
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
    ::std::optional<StreamPreparerFunT> _stream_preparer;
    ::std::optional<StreamSustainerFunT> _stream_sustainer;
    ::std::optional<StreamValuePutterFunT> _stream_value_putter;
    ::std::optional<StreamSeparerFunT> _stream_separer;
    ::std::optional<StreamFinisherFunT> _stream_finisher;
    ::std::optional<bool> _format_independence_flag;
};

// ***********************************************************************

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_preparer(StreamPreparerFunT _) {
    _stream_preparer = _;
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_sustainer(StreamSustainerFunT _) {
    _stream_sustainer = _;
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_value_putter(StreamValuePutterFunT _) {
    _stream_value_putter = _;
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_separer(StreamSeparerFunT _) {
    _stream_separer = _;
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_stream_finisher(StreamFinisherFunT _) {
    _stream_finisher = _;
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_format_independence_flag(bool _) {
    _format_independence_flag = _;
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_null_preparer() {
    _stream_preparer = [](::std::ostream&) {};
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_char_preparer(char c) {
    _stream_preparer = [c](::std::ostream& s) { s << c; };
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_string_preparer(::std::string str) {
    _stream_preparer = [str](::std::ostream& s) { s << str; };
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_null_sustainer() {
    _stream_sustainer = [](::std::ostream&, size_t) {};
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_null_separer() {
    _stream_separer = [](::std::ostream&) {};
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_char_separer(char c) {
    _stream_separer = [c](::std::ostream& s) { s << c; };
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_string_separer(::std::string str) {
    _stream_separer = [str](::std::ostream& s) { s << str; };
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_null_finisher() {
    _stream_finisher = [](::std::ostream&) {};
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_char_finisher(char c) {
    _stream_finisher = [c](::std::ostream& s) { s << c; };
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::set_string_finisher(::std::string str) {
    _stream_finisher = [str](::std::ostream& s) { s << str; };
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_bracket_round() {
    set_char_preparer('(');
    set_char_finisher(')');
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_bracket_curly() {
    set_char_preparer('{');
    set_char_finisher('}');
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_bracket_square() {
    set_char_preparer('[');
    set_char_finisher(']');
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_bracket_angle() {
    set_char_preparer('<');
    set_char_finisher('>');
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_quotation_single() {
    set_char_preparer('\'');
    set_char_finisher('\'');
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_quotation_double() {
    set_char_preparer('"');
    set_char_finisher('"');
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_quotation_back() {
    set_char_preparer('`');
    set_char_finisher('`');
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::in_null() {
    set_char_preparer(' ');
    set_char_finisher(' ');
    set_null_sustainer();
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::like_python_tuple() {
    set_char_preparer('[');
    set_char_finisher(']');
    set_string_separer(", ");
    set_null_sustainer();
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::like_python_list() {
    set_char_preparer('[');
    set_char_finisher(']');
    set_string_separer(", ");
    set_null_sustainer();
    return *this;
}

template <typename _T>
RangeStreamerSettings<_T>& RangeStreamerSettings<_T>::like_python_set() {
    set_char_preparer('{');
    set_char_finisher('}');
    set_string_separer(", ");
    set_null_sustainer();
    return *this;
}

}  // namespace extension::boost

// #######################################################################
// ##  stream_range_impl                                                ##
// #######################################################################

namespace extension::boost {

template <typename RangeT>
void stream_range_impl(
    const ::std::remove_reference_t<RangeT>& range,
    ::std::ostream& os,
    ::std::function<void(::std::ostream&)> stream_preparer,
    ::std::function<void(::std::ostream&, size_t)> stream_sustainer,
    ::std::function<void(::std::ostream&, typename ::boost::range_value<RangeT>::type)> stream_value_putter,
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

}  // namespace extension::boost

// #######################################################################
// ##  RangeStreamer                                                    ##
// #######################################################################

namespace extension::boost {

template <class _RangeT>
class RangeStreamer;

template <class RangeT>
RangeStreamer<RangeT> make_range_streamer(
    RangeT&& range,
    const RangeStreamerSettings<typename ::boost::range_value<RangeT>::type> range_streamer_settings);

template <class _RangeT>
class RangeStreamer {
    static_assert(!::std::is_rvalue_reference_v<_RangeT>,
                  "R must not be a rvalue reference.");
    static_assert(::std::is_lvalue_reference_v<_RangeT> || (!::std::is_reference_v<_RangeT> && !::std::is_const_v<_RangeT>),
                  "R must be of form: `T`, `T&` or `const T&`.");

   public:  // Helper types:
    using RangeT = _RangeT;
    using ValueT = typename ::boost::range_value<RangeT>::type;

   public:  // Factory function:
    friend RangeStreamer<RangeT> make_range_streamer<RangeT>(RangeT&& range, const RangeStreamerSettings<ValueT> range_streamer_settings);

   public:  // API:
    ::std::ostream& stream(::std::ostream& os) const;
    ::std::string str() const;

   private:
    RangeStreamer(::std::add_rvalue_reference_t<RangeT> range, const RangeStreamerSettings<ValueT> range_streamer_settings);

   private:
    const RangeT _range;
    const RangeStreamerSettings<typename ::boost::range_value<RangeT>::type> _range_streamer_settings;
};

template <class _RangeT>
RangeStreamer<_RangeT>::RangeStreamer(
    ::std::add_rvalue_reference_t<RangeT> range,
    const RangeStreamerSettings<ValueT> range_streamer_settings)
    : _range(::std::forward<RangeT>(range)),
      _range_streamer_settings(range_streamer_settings) {
}

template <class _RangeT>
::std::ostream&
RangeStreamer<_RangeT>::stream(::std::ostream& os) const {
    // Defaults:
    using RangeStreamerSettingsT = RangeStreamerSettings<ValueT>;
    const typename RangeStreamerSettingsT::StreamPreparerFunT default_stream_preparer = [](::std::ostream& s) { s << "{"; };
    const typename RangeStreamerSettingsT::StreamSustainerFunT default_stream_sustainer = [](::std::ostream& s, size_t i) { s << i << ":"; };
    const typename RangeStreamerSettingsT::StreamValuePutterFunT default_stream_value_putter = [](::std::ostream& s, ValueT t) { s << t; };
    const typename RangeStreamerSettingsT::StreamSeparerFunT default_stream_separer = [](::std::ostream& s) { s << ","; };
    const typename RangeStreamerSettingsT::StreamFinisherFunT default_stream_finisher = [](::std::ostream& s) { s << "}"; };
    const bool default_format_independence_flag = true;
    // Apply overrules:
    const auto stream_preparer = (_range_streamer_settings._stream_preparer ? *_range_streamer_settings._stream_preparer : default_stream_preparer);
    const auto stream_sustainer = (_range_streamer_settings._stream_sustainer ? *_range_streamer_settings._stream_sustainer : default_stream_sustainer);
    const auto stream_value_putter = (_range_streamer_settings._stream_value_putter ? *_range_streamer_settings._stream_value_putter : default_stream_value_putter);
    const auto stream_separer = (_range_streamer_settings._stream_separer ? *_range_streamer_settings._stream_separer : default_stream_separer);
    const auto stream_finisher = (_range_streamer_settings._stream_finisher ? *_range_streamer_settings._stream_finisher : default_stream_finisher);
    const auto format_independence_flag = (_range_streamer_settings._format_independence_flag ? *_range_streamer_settings._format_independence_flag : default_format_independence_flag);
    // Stream:
    stream_range_impl<RangeT>(
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

template <class _RangeT>
::std::string RangeStreamer<_RangeT>::str() const {
    ::std::ostringstream oss;
    stream(oss);
    return oss.str();
}

template <class RangeT>
RangeStreamer<RangeT> make_range_streamer(RangeT&& range, const RangeStreamerSettings<typename ::boost::range_value<RangeT>::type> range_streamer_settings) {
    return RangeStreamer<RangeT>(::std::forward<RangeT>(range), range_streamer_settings);
}

}  // namespace extension::boost

// #######################################################################
// ##  operator|, Stringifier                                           ##
// #######################################################################

namespace extension::boost {
namespace stream_pragma {

template <class RangeT>
RangeStreamer<RangeT> operator|(RangeT&& range, const RangeStreamerSettings<typename ::boost::range_value<RangeT>::type> range_streamer_settings) {
    return make_range_streamer(::std::forward<RangeT>(range), range_streamer_settings);
}

template <class RangeT>
::std::ostream& operator<<(::std::ostream& os, const RangeStreamer<RangeT>& range_streamer) {
    return range_streamer.stream(os);
}

struct Stringifier {
};

template <class RangeT>
::std::string operator*(const RangeStreamer<RangeT>& range_streamer, Stringifier) {
    return range_streamer.str();
}

template <typename T>
using RSS = RangeStreamerSettings<T>;

}  // namespace stream_pragma
}  // namespace extension::boost
