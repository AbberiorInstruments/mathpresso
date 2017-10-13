// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_MATHPRESSO_P_H
#define _MATHPRESSO_MATHPRESSO_P_H

// [Dependencies]
#include <mathpresso/mathpresso.h>
#include <asmjit/asmjit.h>

#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <memory>

#include <new>

#if defined(_MSC_VER)
#include <windows.h>
# if _MSC_VER >= 1400
#  include <intrin.h>
# endif // _MSC_VER >= 1400
#endif // _MSC_VER

// ============================================================================
// [mathpresso::Architecture]
// ============================================================================

#if (defined(_M_X64  ) || defined(__x86_64) || defined(__x86_64__) || \
     defined(_M_AMD64) || defined(__amd64 ) || defined(__amd64__ ))
# define MATHPRESSO_ARCH_X64          (1)
#else
# define MATHPRESSO_ARCH_X64          (0)
#endif
#if (defined(_M_IX86 ) || defined(__X86__ ) || defined(__i386  ) || \
     defined(__IA32__) || defined(__I86__ ) || defined(__i386__) || \
     defined(__i486__) || defined(__i586__) || defined(__i686__))
# define MATHPRESSO_ARCH_X86          (!MATHPRESSO_ARCH_X64)
#else
# define MATHPRESSO_ARCH_X86          (0)
#endif

#if (defined(_M_ARM  ) || defined(__arm__ ) || defined(__arm) || \
     defined(_M_ARMT ) || defined(__thumb__))
# define MATHPRESSO_ARCH_ARM          (1)
# define MATHPRESSO_ARCH_ARM64        (0)
#else
# define MATHPRESSO_ARCH_ARM          (0)
# define MATHPRESSO_ARCH_ARM64        (0)
#endif

#if MATHPRESSO_ARCH_X86 || MATHPRESSO_ARCH_X64
# define MATHPRESSO_ARCH_64BIT        (MATHPRESSO_ARCH_X64)
# define MATHPRESSO_ARCH_BE           (0)
# define MATHPRESSO_ARCH_LE           (1)
#endif

#if MATHPRESSO_ARCH_ARM || MATHPRESSO_ARCH_ARM64
# define MATHPRESSO_ARCH_64BIT        (MATHPRESSO_ARCH_ARM64)
# define MATHPRESSO_ARCH_BE           (0)
# define MATHPRESSO_ARCH_LE           (1)
#endif

// ============================================================================
// [mathpresso::StdInt]
// ============================================================================

#if defined(__MINGW32__) || defined(__MINGW64__)
# include <sys/types.h>
#endif
#if defined(_MSC_VER) && (_MSC_VER < 1600)
# include <limits.h>
# if !defined(MATHPRESSO_SUPPRESS_STD_TYPES)
#  if (_MSC_VER < 1300)
typedef signed char      int8_t;
typedef signed short     int16_t;
typedef signed int       int32_t;
typedef signed __int64   int64_t;
typedef unsigned char    uint8_t;
typedef unsigned short   uint16_t;
typedef unsigned int     uint32_t;
typedef unsigned __int64 uint64_t;
#  else
typedef __int8           int8_t;
typedef __int16          int16_t;
typedef __int32          int32_t;
typedef __int64          int64_t;
typedef unsigned __int8  uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
#  endif
# endif
# define MATHPRESSO_INT64_C(x) (x##i64)
# define MATHPRESSO_UINT64_C(x) (x##ui64)
#else
# include <stdint.h>
# include <limits.h>
# define MATHPRESSO_INT64_C(x) (x##ll)
# define MATHPRESSO_UINT64_C(x) (x##ull)
#endif

// ============================================================================
// [MathPresso::Likely / Unlikely]
// ============================================================================

#if defined(__GNUC__) || defined(__clang__)
# define MATHPRESSO_LIKELY(exp) __builtin_expect(!!(exp), 1)
# define MATHPRESSO_UNLIKELY(exp) __builtin_expect(!!(exp), 0)
#else
# define MATHPRESSO_LIKELY(exp) exp
# define MATHPRESSO_UNLIKELY(exp) exp
#endif

// ============================================================================
// [MathPresso::Error Handling]
// ============================================================================

//! \internal
#define MATHPRESSO_ASSERT(exp) do { \
  if (!(exp)) \
    ::mathpresso::mpAssertionFailure(__FILE__, __LINE__, #exp); \
  } while (0)

//! \internal
#define MATHPRESSO_ASSERT_NOT_REACHED() do { \
    ::mathpresso::mpAssertionFailure(__FILE__, __LINE__, "Not reached"); \
  } while (0)

//! \internal
#define MATHPRESSO_PROPAGATE(...) \
  do { \
    ::mathpresso::Error _errorValue = __VA_ARGS__; \
    if (MATHPRESSO_UNLIKELY(_errorValue != ::mathpresso::ErrorCode::kErrorOk)) \
      return _errorValue; \
  } while (0)

//! \internal
#define MATHPRESSO_PROPAGATE_(exp, cleanup) \
  do { \
    ::mathpresso::Error _errorCode = (exp); \
    if (MATHPRESSO_UNLIKELY(_errorCode != ::mathpresso::ErrorCode::kErrorOk)) { \
      cleanup \
      return _errorCode; \
    } \
  } while (0)

//! \internal
#define MATHPRESSO_NULLCHECK(ptr) \
  do { \
    if (MATHPRESSO_UNLIKELY(!(ptr))) \
      return MATHPRESSO_TRACE_ERROR(::mathpresso::ErrorCode::kErrorNoMemory); \
  } while (0)

//! \internal
#define MATHPRESSO_NULLCHECK_(ptr, cleanup) \
  do { \
    if (MATHPRESSO_UNLIKELY(!(ptr))) { \
      cleanup \
      return MATHPRESSO_TRACE_ERROR(::mathpresso::ErrorCode::kErrorNoMemory); \
    } \
  } while (0)

#define MATHPRESSO_TRACE_ERROR(error) \
  ::mathpresso::mpTraceError(error)

// ============================================================================
// [mathpresso::]
// ============================================================================

namespace mathpresso
{

	// ============================================================================
	// [mathpresso::Symbols]
	// ============================================================================

	//! Holds MpOperation-objects and variables, and gives an easy way to finding them.
	class Symbols
	{
		using op_ptr_type = std::shared_ptr<MpOperation>;
		using op_map_type = std::map<std::string, std::vector<op_ptr_type>>;

		using var_ptr_type = std::shared_ptr<AstSymbol>;
		using var_map_type = std::map<std::string, var_ptr_type>;

	public:
		std::string name(const std::shared_ptr<MpOperation>  ptr) const;
		op_ptr_type findFunction(const std::string & name, size_t nargs) const;

		//! looks for a MpOperation-Object, where the parameters are complex or real.
		//! if no direct match is found (ie there is no Operation with real parameters),
		//! a conversion from real to complex is return.
		//! generally the first direct match is return, and after that the first match with conversions.
		//! returns nullptr, if no match is found.
		op_ptr_type findFunction(const std::string & name, size_t nargs, bool paramsAreComplex) const;

		std::vector<op_ptr_type> findFunction(const std::string & name) const;
		bool existsFunction(const std::string & name) const;
		var_ptr_type findVariable(const std::string & name) const;

		//! Makes sure that functions with the same name have to have the 
		//! precedence and association.
		void add(const std::string & name, op_ptr_type obj);

		//! add a variable.
		void add(const std::string & name, var_ptr_type obj);
		void remove(const std::string & name);
		std::vector<std::string> names() const;

		void clear();

		std::vector<std::shared_ptr<AstSymbol>> getVariables();
		op_map_type getFunctions() { return _operations; }
	private:
		op_map_type _operations;
		var_map_type _variables;

	};

	// ============================================================================
	// [mathpresso::mpAssertionFailure]
	// ============================================================================

	//! \internal
	//!
	//! MathPresso assertion handler.
	void mpAssertionFailure(const char* file, int line, const char* msg);

	// ============================================================================
	// [mathpresso::mpTraceError]
	// ============================================================================

	Error mpTraceError(Error error);

	// ============================================================================
	// [mathpresso::addBuiltinObjects]
	// ============================================================================
	uint32_t addBuiltinMpObjects(Context * ctx);

	// ============================================================================
	// [mpsl::ErrorReporter]
	// ============================================================================

	//! Error reporter.
	struct ErrorReporter
	{
		//! Set if `OutputLog` is present. MATHPRESSO then checks only this flag to use it.
		static constexpr uint32_t kInternalOptionLog = 0x00010000;

		ErrorReporter(const std::string & body, uint32_t options, OutputLog* log)
			: _body(body.c_str()),
			_len(body.length()),
			_options(options),
			_log(log)
		{
			// These should be handled by MATHPRESSO before the `ErrorReporter` is created.
			MATHPRESSO_ASSERT((log == nullptr && (_options & kInternalOptionLog) == 0) ||
				(log != nullptr && (_options & kInternalOptionLog) != 0));
		}

		// --------------------------------------------------------------------------
		// [Error Handling]
		// --------------------------------------------------------------------------

		bool reportsErrors() const { return (_options & kInternalOptionLog) != 0; }
		bool reportsWarnings() const { return (_options & kOptionVerbose) != 0; }

		void getLineAndColumn(uint32_t position, uint32_t& line, uint32_t& column);

		void onWarning(uint32_t position, const char* fmt, ...);
		void onWarning(uint32_t position, const std::string& msg);

		Error onError(Error error, uint32_t position, const char* fmt, ...);
		Error onError(Error error, uint32_t position, const std::string& msg);

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------
	private:
		const char* _body;
		size_t _len;

		uint32_t _options;
		OutputLog* _log;
	};

	namespace resolver
	{
		using ContextPtr = std::shared_ptr<Context>;
		
		//! find the correct function.
		std::shared_ptr<MpOperation> resolveFunction(ContextPtr ctx, const std::string & name, size_t numargs, bool takesComplex);

		//! find one function with the given name and number of arguments.
		std::shared_ptr<MpOperation> resolveFunction(ContextPtr ctx, const std::string & name, size_t numargs);

		bool existsFunction(ContextPtr ctx, const std::string & nmae);

		//! Get the Variable with the given name.
		std::shared_ptr<AstSymbol> resolveVariable(ContextPtr ctx, const std::string & name, ContextPtr * ctxOut = nullptr);

		std::vector<std::string> separateName(std::string name);


		//! this is a plain version, we only check the Parents for the operation.
		std::string getFunctionName(ContextPtr ctx, std::shared_ptr<MpOperation> operation);

	} // resolver namespace

} // mathpresso namespace

#endif // _MATHPRESSO_MATHPRESSO_P_H
