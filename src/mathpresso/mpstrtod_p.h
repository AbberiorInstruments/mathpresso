// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MPSTRTOD_P_H
#define _MATHPRESSO_MPSTRTOD_P_H

// [Dependencies]
#include <mathpresso/mathpresso_p.h>

// `mathpresso_p.h` includes asmjit, so we can use `ASMJIT_OS_...`.
#if ASMJIT_OS_WINDOWS
# define MATHPRESSO_STRTOD_MSLOCALE
# include <locale.h>
#else
# define MATHPRESSO_STRTOD_XLOCALE
# include <locale.h>
#endif

namespace mathpresso
{

	// ============================================================================
	// [mathpresso::Parser]
	// ============================================================================

	//! Convert String to a double depending on the locale representation.
	struct StrToD
	{
		MATHPRESSO_NO_COPY(StrToD);

#if defined(MATHPRESSO_STRTOD_MSLOCALE)
		StrToD()
		{
			handle = _create_locale(LC_ALL, "C");
		}
		~StrToD()
		{
			_free_locale(handle);
		}

		bool isOk() const
		{
			return handle != nullptr;
		}
		double conv(const char* s, char** end) const
		{
			return _strtod_l(s, end, handle);
		}
	private:
		_locale_t handle;
#elif defined(MATHPRESSO_STRTOD_XLOCALE)
		StrToD()
		{
			handle = newlocale(LC_ALL_MASK, "C", nullptr);
		}
		~StrToD() { freelocale(handle); }

		bool isOk() const
		{
			return handle != nullptr;
		}
		double conv(const char* s, char** end) const
		{
			return strtod_l(s, end, handle);
		}
	private:
		locale_t handle;
#else
		// TODO: this is definitely unsafe.
		bool isOk() const
		{
			return true;
		}
		double conv(const char* s, char** end) const
		{
			return strtod(s, end);
		}
#endif
	};

} // mathpresso namespace

// [Guard]
#endif // _MATHPRESSO_MPSTRTOD_P_H
