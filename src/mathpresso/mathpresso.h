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

// uncomment, if you want to use the same syntax for function-calls.
//#define _REALREWORK

namespace mathpresso
{

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
type(const type& other) = delete; \
type& operator=(const type& other) = delete; \

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
	struct AstSymbol;

	// ============================================================================
	// [mathpresso::TypeDefs]
	// ============================================================================

	//! MathPresso result type (signed integer).
	typedef unsigned int Error;

	//! Prototype of the compiled function generated by MathPresso.
	typedef void(*CompiledFunc)(double* result, void* data);


#ifdef _REALREWORK
	typedef double(*mpFuncDtoD)(double*);
#else
	typedef double(*Arg0Func)();
	typedef double(*Arg1Func)(double);
	typedef double(*Arg2Func)(double, double);
	typedef double(*Arg3Func)(double, double, double);
	typedef double(*Arg4Func)(double, double, double, double);
	typedef double(*Arg5Func)(double, double, double, double, double);
	typedef double(*Arg6Func)(double, double, double, double, double, double);
	typedef double(*Arg7Func)(double, double, double, double, double, double, double);
	typedef double(*Arg8Func)(double, double, double, double, double, double, double, double);
#endif // _REALREWORK

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


	typedef std::complex<double>(*mpFuncpCtoC)(std::complex<double>*);
	typedef double(*mpFuncpCtoD)(std::complex<double>*);
	typedef std::complex<double>(*mpFuncpDtoC)(double*);


	// ============================================================================
	// [mathpresso::ErrorCode]
	// ============================================================================

	//! MathPresso error codes.
	enum ErrorCode
	{
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
	enum Options
	{
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
	enum VariableFlags
	{
		kVariableRW = 0x00000000,
		kVariableRO = 0x00000001,
		kVariableCplx = 0x0000002,
	};

	// ============================================================================
	// [mathpresso::ContextImpl]
	// ============================================================================

	struct ContextImpl
	{
		//! Reference count (atomic).
		uintptr_t _refCount;
	};

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

		op_ptr_type find(const std::string & name, size_t nargs) const;

		//! looks for a MpOperation-Object, where the parameters are complex or real.
		//! if no direct match is found (ie there is no Operation with real parameters),
		//! a conversion from real to complex is return.
		//! generally the first direct match is return, and after that the first match with conversions.
		//! returns nullptr, if no match is found.
		op_ptr_type find(const std::string & name, size_t nargs, bool paramsAreComplex) const;

		std::vector<op_ptr_type> find(const std::string &name) const;

		//! Makes sure that functions with the same name have to have the 
		//! precedence and association.
		void add(const std::string &name, op_ptr_type obj);

		//! add a variable.
		void add(const std::string &name, var_ptr_type obj);

		void remove(const std::string &name);

		std::vector<std::string> names() const;

	private:
		op_map_type _operations;
		var_map_type _variables;

	};

	// ============================================================================
	// [mathpresso::Context]
	// ============================================================================

	struct AstSymbol;
	enum class AstSymbolType;

	//! MathPresso context.
	//!
	//! Context is an environment where you can add/remove constants, variables and
	//! functions. Context is a reference-counted class that is using copy-on-write.
	//! Working with context is reentrant and making weak or deep copy is thread-safe
	//! (reference counting is atomic). It is possible to create one master context
	//! and use it from different threads to compile many expressions.
	//!
	//! reading: thread-safe; writing: not thread-safe!
	struct MATHPRESSO_API Context : public std::enable_shared_from_this<Context>
	{
		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		//! Create a new `Context` instance.
		Context();
		//! Create a new `Context` based on `other`.
		Context(const Context& other);
		//! Destroy the `Context` instance.
		~Context();

		// --------------------------------------------------------------------------
		// [Copy / Reset]
		// --------------------------------------------------------------------------

		//! Delete all symbols.
		Error reset();
		//! Assignment operator.
		Context& operator=(const Context& other);

		// --------------------------------------------------------------------------
		// [Interface]
		// --------------------------------------------------------------------------

		//! Add built-in intrinsics and constants.
		Error addBuiltIns(void);

		//! Add constant to this context.
		Error addConstant(const std::string &name, double value);
		Error addConstant(const std::string &name, std::complex<double> value);
		//! Add variable to this context.
		Error addVariable(const std::string &name, int offset, unsigned int flags = VariableFlags::kVariableRW);

		//! Adding Symbols to the Context, which can contain function calls. See mpoeration.h for more information.
		Error addObject(const std::string &name, std::shared_ptr<MpOperation> obj);

		//! Internal implementation
		Error addSymbol(AstSymbol* &sym, const std::string &name, AstSymbolType type);
		//! Delete symbol from this context.
		Error delSymbol(const std::string &name);

		//! Retrieve a list of all available symbols (Functions, operators and constants)
		Error listSymbols(std::vector<std::string> &syms);
		
		// getter/setter for the subcontexts and parents.
		Error setParent(std::shared_ptr<Context> ctx);
		Error addChild(const std::string & name, std::shared_ptr<Context> ctx);

		std::shared_ptr<Context> getParent()
		{
			return _parent.lock();
		}
		std::shared_ptr<Context> getChild(const std::string & name)
		{
			return _children.at(name);
		}
		
		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		//! Private data not available to the MathPresso public API.
		ContextImpl* _d;

		Symbols _ops;

	protected:
		std::weak_ptr<Context> _parent;
		std::map<std::string, std::shared_ptr<Context>> _children;

	};

	// ============================================================================
	// [mathpresso::Expression]
	// ============================================================================

	//! MathPresso expression.
	struct MATHPRESSO_API Expression
	{
		MATHPRESSO_NO_COPY(Expression);

		// --------------------------------------------------------------------------
		// [Construction / Destruction]
		// --------------------------------------------------------------------------

		//! Create a new `Expression` instance.
		Expression();
		//! Destroy the `Expression` instance.
		~Expression();

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
		Error compile(std::shared_ptr<Context> ctx, const std::string & body, unsigned int options, OutputLog* log = nullptr);

		//! Get whether the `Expression` contains a valid compiled expression.
		bool isCompiled() const;

		//! Get whether the  `expression` can returns a Complex result
		bool isComplex() const { return _isComplex; }

		//! Reset the expression.
		void reset();

		//! Evaluate expression with variable substitutions.
		//!
		//! Returns the result of the evaluated expression, NaN otherwise.
		double evaluate(void* data) const
		{
			double result[2];
			_func(result, data);

			return result[0];
		}

		//! Evaluates expression with variable substitutions.
		//!
		//! This function cannot cope with complex variables, if they are not aligned to
		//! 16 byte boundaries. Use 'alignas(16)' to force the alignment.
		std::complex<double> evaluateComplex(void* data) const
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
	struct MATHPRESSO_API OutputLog
	{
		//! Output type.
		//!
		//! Specifies how much information to return after a program is parsed/compiled.
		enum Message
		{
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
