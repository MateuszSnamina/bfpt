#pragma once

namespace kstate_impl {

// Helper tag classes:
struct CtrFromRange {};
struct CtrFromBuffer {};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static const CtrFromRange ctr_from_range{};
static const CtrFromBuffer ctr_from_buffer{};
#pragma GCC diagnostic pop

} // end of namespace kstate_impl
