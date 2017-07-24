// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPEVAL_P_H
#define _MATHPRESSO_MPEVAL_P_H

// [Dependencies]
#include "./mathpresso_p.h"
#include "./mathpresso/mpcompiler_p.h"
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

	// The `a != a` condition is used to handle NaN values properly. If one of `a`
	// and `b` is NaN the result should be NaN. When `T` isn't a floating point the
	// condition should be removed by C++ compiler.
	template<typename T> MATHPRESSO_INLINE T mpMin(T a, T b) { return (a != a) | (a < b) ? a : b; }
	template<typename T> MATHPRESSO_INLINE T mpMax(T a, T b) { return (a != a) | (a > b) ? a : b; }

	static MATHPRESSO_INLINE double mpGetNan() { static const DoubleBits value = { MATHPRESSO_UINT64_C(0x7FF8000000000000) }; return value.d; }
	static MATHPRESSO_INLINE double mpGetInf() { static const DoubleBits value = { MATHPRESSO_UINT64_C(0x7FF0000000000000) }; return value.d; }
	static MATHPRESSO_INLINE double mpIsNan(double x) { return DoubleBits::fromDouble(x).isNan(); }
	static MATHPRESSO_INLINE bool mpIsInf(double x) { return DoubleBits::fromDouble(x).isInf(); }
	static MATHPRESSO_INLINE bool mpIsFinite(double x) { return DoubleBits::fromDouble(x).isFinite(); }

	static MATHPRESSO_INLINE double mpRound(double x) {
		double y = ::floor(x);
		return y + (x - y >= 0.5 ? double(1.0) : double(0.0));
	}

	static MATHPRESSO_INLINE double mpRoundEven(double x) { return ::rint(x); }
	static MATHPRESSO_INLINE double mpTrunc(double x) { return ::trunc(x); }
	static MATHPRESSO_INLINE double mpFloor(double x) { return ::floor(x); }
	static MATHPRESSO_INLINE double mpCeil(double x) { return ::ceil(x); }

	static MATHPRESSO_INLINE double mpSignBit(double x) { return DoubleBits::fromDouble(x).hi >= 0x80000000U; }
	static MATHPRESSO_INLINE double mpCopySign(double x, double y) {
		DoubleBits bits = DoubleBits::fromDouble(x);
		bits.hi &= 0x7FFFFFFFU;
		bits.hi |= DoubleBits::fromDouble(y).hi & 0x80000000U;
		return bits.d;
	}
	

	static MATHPRESSO_INLINE double mpAvg(double x, double y) { return (x + y) * 0.5; }
	static MATHPRESSO_INLINE double mpMod(double x, double y) { return fmod(x, y); }
	static MATHPRESSO_INLINE double mpAbs(double x) { return ::fabs(x); }
	static MATHPRESSO_INLINE double mpExp(double x) { return ::exp(x); }
	static MATHPRESSO_INLINE double mpPow(double x, double y) { return ::pow(x, y); }

	static MATHPRESSO_INLINE double mpLog(double x) { return ::log(x); }
	static MATHPRESSO_INLINE double mpLog2(double x) { return ::log2(x); }
	static MATHPRESSO_INLINE double mpLog10(double x) { return ::log10(x); }

	static MATHPRESSO_INLINE double mpSqrt(double x) { return ::sqrt(x); }
	static MATHPRESSO_INLINE double mpFrac(double x) { return x - mpFloor(x); }
	static MATHPRESSO_INLINE double mpRecip(double x) { return 1.0 / x; }

	static MATHPRESSO_INLINE double mpSin(double x) { return ::sin(x); }
	static MATHPRESSO_INLINE double mpCos(double x) { return ::cos(x); }
	static MATHPRESSO_INLINE double mpTan(double x) { return ::tan(x); }

	static MATHPRESSO_INLINE double mpSinh(double x) { return ::sinh(x); }
	static MATHPRESSO_INLINE double mpCosh(double x) { return ::cosh(x); }
	static MATHPRESSO_INLINE double mpTanh(double x) { return ::tanh(x); }

	static MATHPRESSO_INLINE double mpAsin(double x) { return ::asin(x); }
	static MATHPRESSO_INLINE double mpAcos(double x) { return ::acos(x); }
	static MATHPRESSO_INLINE double mpAtan(double x) { return ::atan(x); }

	static MATHPRESSO_INLINE double mpAtan2(double y, double x) { return ::atan2(y, x); }
	static MATHPRESSO_INLINE double mpHypot(double x, double y) { return ::hypot(x, y); }

	// Complex functions

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
	
	static MATHPRESSO_INLINE double mpGetReal(std::complex<double> *arg) { return arg->real(); }
	static MATHPRESSO_INLINE double mpGetImag(std::complex<double> *arg) { return arg->imag(); }

	
	static MATHPRESSO_INLINE std::complex<double> mpLog2C(std::complex<double> *x) { return std::log(*x) / std::log(2); }
	
	static MATHPRESSO_INLINE std::complex<double> mpRecipC(std::complex<double> *x) { return 1.0 / *x; }

	static MATHPRESSO_INLINE double mpAbsC(std::complex<double> *x) { return std::abs(*x); }

	static MATHPRESSO_INLINE std::complex<double> mpAvgC(std::complex<double> *x) { return (x[0] + x[1])*0.5; }

	// Functions for the Optimizer
	static double mpAddOptD(double l, double r) { return l+ r; }
	static std::complex<double> mpAddOptC(std::complex<double> *r) { return r[0] + r[1]; }

	static double mpSubOptD(double l, double r) { return l - r; }
	static std::complex<double> mpSubOptC(std::complex<double> *r) { return r[0] - r[1]; }

	static double mpMulOptD(double l, double r) { return l * r; }
	static std::complex<double> mpMulOptC(std::complex<double> *r) { return r[0] * r[1]; }

	static double mpDivOptD(double l, double r) { return l / r; }
	static std::complex<double> mpDivOptC(std::complex<double> *r) { return r[0] / r[1]; }

	//! The second parameter contains the parameters from left to right.
	//! ie: a binary function a-b --> {param[0] - param[1]}
	typedef JitVar (*mpAsmFunc)(JitCompiler*, JitVar*);

	// the emitters for the different operators.
	static JitVar compileAddD(JitCompiler* jc, JitVar* vars ) 
	{
		jc->cc->addsd(vars[0].getXmm(), vars[1].getXmm());
		return vars[0];
	}
	static JitVar compileAddC(JitCompiler* jc, JitVar* vars) 
	{
		jc->cc->addpd(vars[0].getXmm(), vars[1].getXmm());
		return vars[0];
	}

	static JitVar compileSubD(JitCompiler* jc, JitVar* vars) 
	{
		jc->cc->subsd(vars[0].getXmm(), vars[1].getXmm());
		return vars[0];
	}
	static JitVar compileSubC(JitCompiler* jc, JitVar* vars) 
	{
		jc->cc->subpd(vars[0].getXmm(), vars[1].getXmm());
		return vars[0];
	}
	
	static JitVar compileMulD(JitCompiler* jc, JitVar* vars) 
	{
		jc->cc->mulsd(vars[0].getXmm(), vars[1].getXmm());
		return vars[0];
	}
	
	static JitVar compileMulC(JitCompiler* jc, JitVar* vars) 
	{
		JitVar negateImag = jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000));
		JitVar ret(jc->cc->newXmmPd(), JitVar::FLAG_NONE);
		jc->cc->movapd(ret.getXmm(), vars[0].getXmm());
		jc->cc->mulpd(ret.getXmm(), vars[1].getXmm());
		jc->cc->shufpd(vars[1].getXmm(), vars[1].getXmm(), 1);
		jc->cc->pxor(vars[1].getXmm(), negateImag.getMem());
		jc->cc->mulpd(vars[0].getXmm(), vars[1].getXmm());
		jc->cc->hsubpd(ret.getXmm(), vars[0].getXmm());
		return ret;
	}

	static JitVar compileDivD(JitCompiler* jc, JitVar* vars) 
	{
		jc->cc->divsd(vars[0].getXmm(), vars[1].getXmm());
		return vars[0];
	}
	static JitVar compileDivC(JitCompiler* jc, JitVar* vars) 
	{
		JitVar negateImag = jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000));
		JitVar ret(jc->cc->newXmmPd(), JitVar::FLAG_NONE);
		jc->cc->pxor(vars[1].getXmm(), negateImag.getMem());

		jc->cc->movapd(ret.getXmm(), vars[0].getXmm());
		jc->cc->mulpd(ret.getXmm(), vars[1].getXmm());
		jc->cc->shufpd(vars[1].getXmm(), vars[1].getXmm(), 1);
		jc->cc->pxor(vars[1].getXmm(), negateImag.getMem());
		jc->cc->mulpd(vars[0].getXmm(), vars[1].getXmm());
		jc->cc->hsubpd(ret.getXmm(), vars[0].getXmm());

		jc->cc->mulpd(vars[1].getXmm(), vars[1].getXmm());
		jc->cc->haddpd(vars[1].getXmm(), vars[1].getXmm());
		jc->cc->divpd(ret.getXmm(), vars[1].getXmm());
		return ret;
	}

} // mathpresso namespace

#endif // _MATHPRESSO_MPEVAL_P_H
