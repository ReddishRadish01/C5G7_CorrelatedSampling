#pragma once

// Compile-time tally switches.
//
// Pin-power benchmark mode: collision estimator of the fission reaction rate,
// weight * Sigma_f / Sigma_collision. This intentionally excludes nu.
#define C5G7_TALLY_PINWISE_FISSION_RATE

// Alternative source-distribution mode: expected fission-neutron source
// weight, nu * weight * Sigma_f / Sigma_collision.
// #define C5G7_TALLY_FISSION_SOURCE_DIST

#if defined(C5G7_TALLY_PINWISE_FISSION_RATE) && defined(C5G7_TALLY_FISSION_SOURCE_DIST)
#error Select only one tally quantity: pin-wise fission rate or fission-source distribution.
#endif

// When enabled, main.cu writes every assembly and every axial tally cell to a
// separate pinwise_Tally_*.txt file after the final tally fetch.
#define C5G7_DUMP_FULL_CORE_PINWISE_TALLY

#if defined(C5G7_DUMP_FULL_CORE_PINWISE_TALLY) && \
	!defined(C5G7_TALLY_PINWISE_FISSION_RATE) && \
	!defined(C5G7_TALLY_FISSION_SOURCE_DIST)
#error Full-core pin-wise dump requires one tally quantity.
#endif
