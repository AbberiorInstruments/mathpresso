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
    if (MATHPRESSO_UNLIKELY(_errorValue != ::mathpresso::kErrorOk)) \
      return _errorValue; \
  } while (0)

//! \internal
#define MATHPRESSO_PROPAGATE_(exp, cleanup) \
  do { \
    ::mathpresso::Error _errorCode = (exp); \
    if (MATHPRESSO_UNLIKELY(_errorCode != ::mathpresso::kErrorOk)) { \
      cleanup \
      return _errorCode; \
    } \
  } while (0)

//! \internal
#define MATHPRESSO_NULLCHECK(ptr) \
  do { \
    if (MATHPRESSO_UNLIKELY(!(ptr))) \
      return MATHPRESSO_TRACE_ERROR(::mathpresso::kErrorNoMemory); \
  } while (0)

//! \internal
#define MATHPRESSO_NULLCHECK_(ptr, cleanup) \
  do { \
    if (MATHPRESSO_UNLIKELY(!(ptr))) { \
      cleanup \
      return MATHPRESSO_TRACE_ERROR(::mathpresso::kErrorNoMemory); \
    } \
  } while (0)

#define MATHPRESSO_TRACE_ERROR(error) \
  ::mathpresso::mpTraceError(error)

// ============================================================================
// [mathpresso::]
// ============================================================================

namespace mathpresso {

// ============================================================================
// [Reuse]
// ============================================================================

// Reuse these classes - we depend on asmjit anyway and these are internal.
using asmjit::StringBuilder;
using asmjit::StringBuilderTmp;

using asmjit::Zone;
using asmjit::ZoneHeap;
using asmjit::ZoneVector;

// ============================================================================
// [mathpresso::OpFlags]
// ============================================================================

//! Operator flags.
enum OpFlags {
  //! The operator has one parameter (unary node).
  kOpFlagUnary         = 0x00000001,
  //! The operator has two parameters (binary node).
  kOpFlagBinary        = 0x00000002,
  //! The operator is an intrinsic function.
  kOpFlagIntrinsic     = 0x00000004,
  //! The operator has right-to-left associativity.
  kOpFlagRightToLeft   = 0x00000008,

  //! The operator does an assignment to a variable.
  kOpFlagAssign        = 0x00000010,

  //! The operator 3 parameters (ternary node).
  kOpFlagTernary	   = 0x00000020,

  //! this Operator is complex -> complex.
  kOpFlagComplexToComplex = 0x00001000,
  //! this Operator is real -> real
  kOpFlagRealToReal = 0x00002000,
  //! this Operator is complex -> real
  kOpFlagComplexToReal = 0x0004000,
  //! this Operator is real -> complex
  kOpFlagRealToComplex = 0x0008000,

  //! The operator performs an arithmetic operation.
  kOpFlagArithmetic    = 0x00000100,
  //! The operator performs a conditional operation.
  kOpFlagCondition     = 0x00000200,
  //! The operator performs a floating-point rounding.
  kOpFlagRounding      = 0x00000400,
  //! The operator calculates a trigonometric function.
  kOpFlagTrigonometric = 0x00000800,

  // set if there is an object for this function.
  _kOpFlagHasobject = 0x00010000,

  kOpFlagNopIfLZero    = 0x10000000,
  kOpFlagNopIfRZero    = 0x20000000,
  kOpFlagNopIfLOne     = 0x40000000,
  kOpFlagNopIfROne     = 0x80000000,

  kOpFlagNopIfZero     = kOpFlagNopIfLZero | kOpFlagNopIfRZero,
  kOpFlagNopIfOne      = kOpFlagNopIfLOne  | kOpFlagNopIfROne
};

// ============================================================================
// [mpsl::InternalConsts]
// ============================================================================

enum InternalConsts {
  kInvalidSlot = 0xFFFFFFFFU
};

// ============================================================================
// [mpsl::InternalOptions]
// ============================================================================

//! \internal
//!
//! Compilation options MATHPRESSO uses internally.
enum InternalOptions {
  //! Set if `OutputLog` is present. MATHPRESSO then checks only this flag to use it.
  kInternalOptionLog = 0x00010000
};

// ============================================================================
// [mathpresso::mpAssertionFailure]
// ============================================================================

//! \internal
//!
//! MathPresso assertion handler.
MATHPRESSO_NOAPI void mpAssertionFailure(const char* file, int line, const char* msg);

// ============================================================================
// [mathpresso::mpTraceError]
// ============================================================================

MATHPRESSO_NOAPI Error mpTraceError(Error error);


// ============================================================================
// [mathpresso::addBuiltinObjects]
// ============================================================================
struct Context;
uint32_t addBuiltinMpObjects(Context * ctx);

// ============================================================================
// [mathpresso::StringRef]
// ============================================================================

//! String reference (pointer to string data data and length).
//!
//! NOTE: MATHPRESSO always provides NULL terminated string with length known. On the
//! other hand MATHPRESSO doesn't require NULL terminated strings when passed to MATHPRESSO
//! APIs.
struct StringRef {
  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  MATHPRESSO_INLINE StringRef()
    : _data(nullptr),
      _length(0) {}

  explicit MATHPRESSO_INLINE StringRef(const char* data)
    : _data(data),
      _length(::strlen(data)) {}

  MATHPRESSO_INLINE StringRef(const char* data, size_t len)
    : _data(data),
      _length(len) {}

  // --------------------------------------------------------------------------
  // [Reset / Setup]
  // --------------------------------------------------------------------------

  MATHPRESSO_INLINE void reset() {
    set(nullptr, 0);
  }

  MATHPRESSO_INLINE void set(const char* s) {
    set(s, ::strlen(s));
  }

  MATHPRESSO_INLINE void set(const char* s, size_t len) {
    _data = s;
    _length = len;
  }

  // --------------------------------------------------------------------------
  // [Accessors]
  // --------------------------------------------------------------------------

  //! Get the string data.
  MATHPRESSO_INLINE const char* getData() const { return _data; }
  //! Get the string length.
  MATHPRESSO_INLINE size_t getLength() const { return _length; }

  // --------------------------------------------------------------------------
  // [Eq]
  // --------------------------------------------------------------------------

  MATHPRESSO_INLINE bool eq(const char* s) const {
    const char* a = _data;
    const char* aEnd = a + _length;

    while (a != aEnd) {
      if (*a++ != *s++)
        return false;
    }

    return *s == '\0';
  }

  MATHPRESSO_INLINE bool eq(const char* s, size_t len) const {
    return len == _length && ::memcmp(_data, s, len) == 0;
  }

  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  const char* _data;
  size_t _length;
};

// ============================================================================
// [mpsl::ErrorReporter]
// ============================================================================

//! Error reporter.
struct ErrorReporter {
  MATHPRESSO_INLINE ErrorReporter(const char* body, size_t len, uint32_t options, OutputLog* log)
    : _body(body),
      _len(len),
      _options(options),
      _log(log) {

    // These should be handled by MATHPRESSO before the `ErrorReporter` is created.
    MATHPRESSO_ASSERT((log == nullptr && (_options & kInternalOptionLog) == 0) ||
                      (log != nullptr && (_options & kInternalOptionLog) != 0) );
  }

  // --------------------------------------------------------------------------
  // [Error Handling]
  // --------------------------------------------------------------------------

  MATHPRESSO_INLINE bool reportsErrors() const { return (_options & kInternalOptionLog) != 0; }
  MATHPRESSO_INLINE bool reportsWarnings() const { return (_options & kOptionVerbose) != 0; }

  void getLineAndColumn(uint32_t position, uint32_t& line, uint32_t& column);

  void onWarning(uint32_t position, const char* fmt, ...);
  void onWarning(uint32_t position, const StringBuilder& msg);

  Error onError(Error error, uint32_t position, const char* fmt, ...);
  Error onError(Error error, uint32_t position, const StringBuilder& msg);

  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  const char* _body;
  size_t _len;

  uint32_t _options;
  OutputLog* _log;
};

} // mathpresso namespace

#endif // _MATHPRESSO_MATHPRESSO_P_H
