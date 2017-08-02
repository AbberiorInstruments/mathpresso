// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MATHPRESSO_H
#define _MATHPRESSO_H

#include <stdlib.h>
#include <iostream>
#include <complex>
#include <vector>
#include <map>
#include <memory>
#include <algorithm>

#if !defined(_MSC_VER)
#include <stdint.h>
#endif

namespace mathpresso {

// ============================================================================
// [mathpresso::Configuration]
// ============================================================================

// EMBED implies STATIC.
#if defined(MATHPRESSO_EMBED) && !defined(MATHPRESSO_STATIC)
# define MATHPRESSO_STATIC
#endif

// ============================================================================
// [mathpresso::PPDefs]
// ============================================================================

//! \def MATHPRESSO_API
//!
//! Mathpresso API decorator.
#if !defined(MATHPRESSO_API)
# if defined(MATHPRESSO_STATIC)
#  define MATHPRESSO_API
# elif defined(_WINDOWS)
#  if defined(__GNUC__) || defined(__clang__) && !defined(__MINGW32__)
#   if defined(MATHPRESSO_EXPORTS)
#    define MATHPRESSO_API __attribute__((__dllexport__))
#   else
#    define MATHPRESSO_API __attribute__((__dllimport__))
#   endif
#  else
#   if defined(MATHPRESSO_EXPORTS)
#    define MATHPRESSO_API __declspec(dllexport)
#   else
#    define MATHPRESSO_API __declspec(dllimport)
#   endif
#  endif
# else
#  if defined(__clang__) || defined(__GNUC__)
#   define MATHPRESSO_API __attribute__((__visibility__("default")))
#  endif
# endif
#endif

//! \def MATHPRESSO_NOAPI
//!
//! Mathpresso hidden API decorator.
#define MATHPRESSO_NOAPI

//! \def MATHPRESSO_INLINE
//!
//! Mathpresso inline decorator.
#if defined(__clang__)
# define MATHPRESSO_INLINE inline __attribute__((__always_inline__, __visibility__("hidden")))
#elif defined(__GNUC__)
# define MATHPRESSO_INLINE inline __attribute__((__always_inline__))
#elif defined(_MSC_VER)
# define MATHPRESSO_INLINE __forceinline
#else
# define MATHPRESSO_INLINE inline
#endif

#define MATHPRESSO_NO_COPY(type) \
private: \
  MATHPRESSO_INLINE type(const type& other); \
  MATHPRESSO_INLINE type& operator=(const type& other); \
public:

//! Get an offset of `field` in a struct `type`.
#define MATHPRESSO_OFFSET(type, field) \
  ((int)(size_t) ((const char*) &((const type*)0x10)->field) - 0x10)

#define MATHPRESSO_ARRAY_SIZE(array) \
  (sizeof(array) / sizeof(array[0]))

// ============================================================================
// [Forward Declarations]
// ============================================================================

struct OutputLog;
struct Expression;
class MpOperation;

// ============================================================================
// [mathpresso::TypeDefs]
// ============================================================================

//! MathPresso result type (signed integer).
typedef unsigned int Error;

//! Prototype of the compiled function generated by MathPresso.
typedef void (*CompiledFunc)(double* result, void* data);

typedef double (*Arg0Func)();
typedef double (*Arg1Func)(double);
typedef double (*Arg2Func)(double, double);
typedef double (*Arg3Func)(double, double, double);
typedef double (*Arg4Func)(double, double, double, double);
typedef double (*Arg5Func)(double, double, double, double, double);
typedef double (*Arg6Func)(double, double, double, double, double, double);
typedef double (*Arg7Func)(double, double, double, double, double, double, double);
typedef double (*Arg8Func)(double, double, double, double, double, double, double, double);

template <std::complex<double>(*T)(const std::complex<double> &)>
std::complex<double> mpFuncCtoC1(std::complex<double> *x)
{
	return T(*x);
}

template <std::complex<double>(*T)(const std::complex<double> &, const std::complex<double> &)>
std::complex<double> mpFuncCtoC2(std::complex<double> *x)
{
	return T(x[0], x[1]);
}

template <std::complex<double>(*T)(const std::complex<double> &, const std::complex<double> &, const std::complex<double> &)>
std::complex<double> mpFuncCtoC3(std::complex<double> *x)
{
	return T(x[0], x[1], x[2]);
}

typedef std::complex<double> (*mpFuncpCtoC)(std::complex<double>*);
typedef double (*mpFuncpCtoD)(std::complex<double>*);
typedef std::complex<double>(*mpFuncpDtoC)(double*);


// ============================================================================
// [mathpresso::ErrorCode]
// ============================================================================

//! MathPresso error codes.
enum ErrorCode {
  //! No error.
  kErrorOk = 0,
  //! No memory.
  kErrorNoMemory = 1,
  //! Invalid argument.
  kErrorInvalidArgument,
  //! Invalid state.
  kErrorInvalidState,

  //! No expression was given.
  kErrorNoExpression,
  //! Invalid syntax.
  kErrorInvalidSyntax,

  //! Symbol not found.
  kErrorSymbolNotFound,
  //! Symbol already exists.
  kErrorSymbolAlreadyExists
};

// ============================================================================
// [mathpresso::Options]
// ============================================================================

//! MathPresso options.
enum Options {
  //! No options.
  kNoOptions = 0x00000000,

  //! Show messages and warnings.
  kOptionVerbose = 0x0001,
  //! Debug AST (shows initial and final AST).
  kOptionDebugAst = 0x0002,
  //! Debug assembly generated.
  kOptionDebugAsm = 0x0008,

  //! Do not use SSE4.1 instruction set even if CPU supports it.
  //!
  //! NOTE: This is used during testing to ensure that all code-paths produce
  //! the same results regardless of the highest instruction set used. Since
  //! SSE4.1 is the most beneficial instruction set for MathPresso there is
  //! only this option (MathPresso doesn't use SSE3 and SSSE3 at the moment).
  kOptionDisableSSE4_1 = 0x4000,

  //! \internal
  //!
  //! Mask of all accessible options, MathPresso uses also \ref InternalOptions
  //! that should not collide with \ref Options.
  _kOptionsMask = 0xFFFF
};

// ============================================================================
// [mathpresso::VariableFlags]
// ============================================================================

//! Variable flags.
enum VariableFlags {
  kVariableRW = 0x00000000,
  kVariableRO = 0x00000001,
  kVariableCplx = 0x0000002,
};

// ============================================================================
// [mathpresso::FunctionFlags]
// ============================================================================

enum FunctionFlags {
  //! Function has 0 arguments.
  kFunctionArg0 = 0x00000000,
  //! Function has 1 argument.
  kFunctionArg1 = 0x00000001,
  //! Function has 2 arguments.
  kFunctionArg2 = 0x00000002,
  //! Function has 3 arguments.
  kFunctionArg3 = 0x00000003,
  //! Function has 4 arguments.
  kFunctionArg4 = 0x00000004,
  //! Function has 5 arguments.
  kFunctionArg5 = 0x00000005,
  //! Function has 6 arguments.
  kFunctionArg6 = 0x00000006,
  //! Function has 7 arguments.
  kFunctionArg7 = 0x00000007,
  //! Function has 8 arguments.
  kFunctionArg8 = 0x00000008,

  //! \internal
  _kFunctionArgMask = 0x0000000F,

  kFunctionReturnsComplex = 0x00000010,
  kFunctionTakesComplex = 0x00000020,

  //! Function cannot be optimized out as the pointer carries state information
  //! that can change
  kFunctionHasState = 0x00000040,

  //! The first argument of the function is the `data` pointer passed to the
  //! evaluate function. This is a hidden parameter that is not accessible
  //! within the expression itself.
  kFunctionFirstArgData = 0x10000000,

  //! Function doesn't have side-effects and can be evaluated (i.e. optimized
  //! out) during a constant folding phase.
  kFunctionNoSideEffects = 0x20000000,

  _kFunctionHasObject = 0x20000000
};

// ============================================================================
// [mathpresso::ContextImpl]
// ============================================================================

struct ContextImpl {
  //! Reference count (atomic).
  uintptr_t _refCount;
};

// ============================================================================
// [mathpresso::Operations]
// ============================================================================

typedef std::map<std::pair<std::string, size_t>, std::shared_ptr<MpOperation>> symbolMap;

//! Holds MpOperation-objects and gives an easy way to finding them.
class Operations
{
public:

	//! Get The name of a function, where you have a pointer to MpOperation.
	std::string findName(MpOperation* ptr) const;

	//! get the mame of a function where have a shared_ptr to MpOperation
	std::string findName(std::shared_ptr<MpOperation> ptr) const;	

	//! Get a shared_ptr to a MpOperation, by its name and number of arguments.
	std::shared_ptr<MpOperation> getOperation(std::string name, size_t numArgs) const;	

	//! Get all MpOperations with a specific name, regardless of the number of arguments.
	std::vector<std::shared_ptr<MpOperation>> findOperations(std::string name) const;	

	//! Add a MpOperation-Object to the internal map.
	void addOperation(std::string name, std::shared_ptr<MpOperation> obj);

	//! Check whether a MpOperation-Object with a known name and number of arguments exists.
	bool hasOperation(std::string name, size_t numArgs) const;

private:
	symbolMap _symbols;
};

// ============================================================================
// [mathpresso::Context]
// ============================================================================

struct AstSymbol;

//! MathPresso context.
//!
//! Context is an environment where you can add/remove constants, variables and
//! functions. Context is a reference-counted class that is using copy-on-write.
//! Working with context is reentrant and making weak or deep copy is thread-safe
//! (reference counting is atomic). It is possible to create one master context
//! and use it from different threads to compile many expressions.
struct Context {
  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  //! Create a new `Context` instance.
  MATHPRESSO_API Context();
  //! Create a new `Context` based on `other`.
  MATHPRESSO_API Context(const Context& other);
  //! Destroy the `Context` instance.
  MATHPRESSO_API ~Context();

  // --------------------------------------------------------------------------
  // [Copy / Reset]
  // --------------------------------------------------------------------------

  //! Delete all symbols.
  MATHPRESSO_API Error reset();
  //! Assignement operator.
  MATHPRESSO_API Context& operator=(const Context& other);

  // --------------------------------------------------------------------------
  // [Interface]
  // --------------------------------------------------------------------------

  //! Add built-in intrinsics and constants.
  MATHPRESSO_API Error addBuiltIns(void);

  //! Add constant to this context.
  MATHPRESSO_API Error addConstant(const char* name, double value);
  MATHPRESSO_API Error addConstant(const char * name, std::complex<double> value);
  //! Add variable to this context.
  MATHPRESSO_API Error addVariable(const char* name, int offset, unsigned int flags = kVariableRW);

  //! Adding Operations to the Context, which can contain function calls. See mpoeration.h for more information.
  MATHPRESSO_API Error addObject(std::string name, std::shared_ptr<MpOperation> obj);

  //! Internal implementation
  MATHPRESSO_API Error addSymbol(AstSymbol* &sym, const char * name, int type);
  //! Delete symbol from this context.
  MATHPRESSO_API Error delSymbol(const char* name);

  //! Retrieve a list of all available symbols
  MATHPRESSO_API Error listSymbols(std::vector<std::string> &syms);
  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  //! Private data not available to the MathPresso public API.
  ContextImpl* _d;



  Operations _ops;
};

// ============================================================================
// [mathpresso::Expression]
// ============================================================================

//! MathPresso expression.
struct Expression {
  MATHPRESSO_NO_COPY(Expression)

  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  //! Create a new `Expression` instance.
  MATHPRESSO_API Expression();
  //! Destroy the `Expression` instance.
  MATHPRESSO_API ~Expression();

  // --------------------------------------------------------------------------
  // [Interface]
  // --------------------------------------------------------------------------

  //! Parse and compile a given expression.
  //!
  //! \param ctx MathPresso's \ref Context to use.
  //! \param body Expression to parse and compile.
  //! \param options MathPresso options (flags), see \ref Options.
  //! \param log Used to catch messages a parser, optimizer, and compiler may
  //!        generate
  //!
  //! Returns MathPresso's error code, see \ref Error.
  MATHPRESSO_API Error compile(const Context& ctx, const char* body, unsigned int options, OutputLog* log = nullptr);

  //! Get whether the `Expression` contains a valid compiled expression.
  MATHPRESSO_API bool isCompiled() const;

  //! Get whether the  `expression` can returns a Complex result
  MATHPRESSO_API bool isComplex() const { return _isComplex; }

  //! Reset the expression.
  MATHPRESSO_API void reset();

  //! Evaluate expression with variable substitutions.
  //!
  //! Returns the result of the evaluated expression, NaN otherwise.
  MATHPRESSO_INLINE double evaluate(void* data) const
  {
	  double result[2];
	  _func(result, data);

	  return result[0];
  }

  //! Evaluates expression with variable substitutions.
  //!
  //! This function cannot cope with complex variables, if they are not aligned to
  //! 16 byte boundaries. Use 'alignas(16)' to force the alignment.
  MATHPRESSO_INLINE std::complex<double> evaluateComplex(void* data) const 
  {
	  double result[2] = { 0, 0 };
	  _func(result, data);

	  return std::complex<double>(result[0], result[1]);
  }


  // --------------------------------------------------------------------------
  // [Members]
  // --------------------------------------------------------------------------

  //! Compiled function.
  CompiledFunc _func;


  //! True, if the result of the function could be a complex number
  bool _isComplex = false;
};


// ============================================================================
// [mpsl::OutputLog]
// ============================================================================

//! Interface that can be used to catch compiler warnings and errors.
struct MATHPRESSO_API OutputLog {
  //! Output type.
  //!
  //! Specifies how much information to return after a program is parsed/compiled.
  enum Message {
    //! Error message.
    kMessageError = 0,
    //! warning.
    kMessageWarning,
    //! AST initial.
    kMessageAstInitial,
    //! AST final.
    kMessageAstFinal,
    //! Machine code.
    kMessageAsm
  };

  // --------------------------------------------------------------------------
  // [Construction / Destruction]
  // --------------------------------------------------------------------------

  OutputLog();
  virtual ~OutputLog();

  // --------------------------------------------------------------------------
  // [Interface]
  // --------------------------------------------------------------------------

  virtual void log(unsigned int type, unsigned int line, unsigned int column, const char* message, size_t len) = 0;
};

} // mathpresso namespace

#endif // _MATHPRESSO_H
