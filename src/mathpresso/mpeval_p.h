// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPEVAL_P_H
#define _MATHPRESSO_MPEVAL_P_H

// [Dependencies]
#include <mathpresso/mathpresso_p.h>
#include <mathpresso/mpcompiler_p.h>
#include <complex>
#include <iostream>

namespace mathpresso {

	// ============================================================================
	// [mathpresso::Evaluation]
	// ============================================================================

	//! DP-FP binary representation and utilities.
	union DoubleBits {
		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		static MATHPRESSO_INLINE DoubleBits fromDouble(double val) { DoubleBits u; u.d = val; return u; }
		static MATHPRESSO_INLINE DoubleBits fromUInt64(uint64_t val) { DoubleBits u; u.u = val; return u; }

		MATHPRESSO_INLINE bool isNan() const { return ((hi & 0x7FF00000U)) == 0x7FF00000U && ((hi & 0x000FFFFFU) | lo) != 0x00000000U; }
		MATHPRESSO_INLINE void setNan() { u = MATHPRESSO_UINT64_C(0x7FF8000000000000); }

		MATHPRESSO_INLINE bool isInf() const { return (hi & 0x7FFFFFFFU) == 0x7FF00000U && lo == 0x00000000U; }
		MATHPRESSO_INLINE void setInf() { u = MATHPRESSO_UINT64_C(0x7FF0000000000000); }

		MATHPRESSO_INLINE bool isFinite() const { return (hi & 0x7FF00000U) != 0x7FF00000U; }

		//! Value as uint64_t.
		uint64_t u;
		//! Value as `double`.
		double d;

#if MATHPRESSO_ARCH_LE
		struct { uint32_t lo; uint32_t hi; };
#else
		struct { uint32_t hi; uint32_t lo; };
#endif
	};

	//! DP-FP binary representation and utilities for complex numbers.
	// If either the real or the imaginary part of the number is NaN or Inf, the whole Complex is NaN or Inf.
	union DoubleBitsComp {
		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		static MATHPRESSO_INLINE DoubleBitsComp fromDouble(double val0, double val1) { DoubleBitsComp u; u.d[0] = val0; u.d[1] = val1; return u; }
		static MATHPRESSO_INLINE DoubleBitsComp fromUInt64(uint64_t val0, uint64_t val1) { DoubleBitsComp u; u.u[0] = val0; u.u[1] = val1; return u; }
		static MATHPRESSO_INLINE DoubleBitsComp fromDoubleComplex(std::complex<double> val) {
			DoubleBitsComp u;
			u.d[0] = val.real();
			u.d[1] = val.imag();
			return u;
		}

		MATHPRESSO_INLINE bool isNan() const {
			return (((hi0 & 0x7FF00000U)) == 0x7FF00000U && ((hi0 & 0x000FFFFFU) | lo0) != 0x00000000U)
				|| (((hi1 & 0x7FF00000U)) == 0x7FF00000U && ((hi1 & 0x000FFFFFU) | lo1) != 0x00000000U);
		}
		MATHPRESSO_INLINE void setNan() { u[0] = u[1] = MATHPRESSO_UINT64_C(0x7FF8000000000000); }

		MATHPRESSO_INLINE bool isInf() const {
			return ((hi0 & 0x7FFFFFFFU) == 0x7FF00000U && lo0 == 0x00000000U)
				|| ((hi1 & 0x7FFFFFFFU) == 0x7FF00000U && lo1 == 0x00000000U);
		}
		MATHPRESSO_INLINE void setInf() { u[0] = u[1] = MATHPRESSO_UINT64_C(0x7FF0000000000000); }

		MATHPRESSO_INLINE bool isFinite() const {
			return (hi0 & 0x7FF00000U) != 0x7FF00000U
				&& (hi1 & 0x7FF00000U) != 0x7FF00000U;
		}

		//! Value as uint64_t.
		uint64_t u[2];
		//! Value as `double`.
		double d[2];

#if MATHPRESSO_ARCH_LE
		struct { uint32_t lo0; uint32_t hi0; uint32_t lo1; uint32_t hi1; };
#else
		struct { uint32_t hi0; uint32_t lo0; uint32_t hi1; uint32_t lo1; };
#endif
	};

	static MATHPRESSO_INLINE double mpGetNan() {  return std::numeric_limits<double>::quiet_NaN(); }
	static MATHPRESSO_INLINE double mpGetInf() {  return std::numeric_limits<double>::infinity(); }
	

	//! Used to call a cpp-function from within the assembler.
	//! the result can be read from data at index 0, arguments are at index 1, 2, ...
	//! Eventually there is a better way, but unless i understand how to pass structs, i cannot provide one.
	static MATHPRESSO_INLINE void mpWrapCtoC(std::complex<double>(*ptr)(std::complex<double> *), std::complex<double>* data) {
		data[0] = ptr(data + 1);
	}

	// function from double to complex
	static MATHPRESSO_INLINE void mpWrapDtoC(std::complex<double>(*ptr)(double *), double* data, std::complex<double>* ret) {
		*ret = ptr(data);
	}

} // mathpresso namespace

#endif // _MATHPRESSO_MPEVAL_P_H
