// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include <mathpresso/mathpresso_p.h>
#include <mathpresso/mpast_p.h>
#include <mathpresso/mpatomic_p.h>
#include <mathpresso/mpcompiler_p.h>
#include <mathpresso/mpeval_p.h>
#include <mathpresso/mpoptimizer_p.h>
#include <mathpresso/mpparser_p.h>
#include <mathpresso/mptokenizer_p.h>
#include <mathpresso/mpoperation.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <map>

namespace mathpresso
{

	// ============================================================================
	// [mathpresso::mpAssertionFailure]
	// ============================================================================

	void mpAssertionFailure(const char* file, int line, const char* msg)
	{
		fprintf(stderr,
				"[mathpresso] Assertion failed at %s (line %d):\n"
				"[mathpresso] %s\n", file, line, msg);

		::abort();
	}

	// ============================================================================
	// [mathpresso::mpTraceError]
	// ============================================================================

	Error mpTraceError(Error error)
	{
		return error;
	}

	// ============================================================================
	// [mathpresso::mpDummyFunc]
	// ============================================================================

	//! \internal
	//!
	//! Used instead of NULL in `Expression::_func`.
	//!
	//! Returns NaN.
	static void mpDummyFunc(double* result, void*)
	{
		result[0] = mpGetNan();
		result[1] = mpGetNan();
	}

	// ============================================================================
	// [mathpresso::Context - Copy / Reset]
	// ============================================================================

	Error Context::reset()
	{
		_parent.reset();
		_children.clear();
		_symbols->clear();
		return ErrorCode::kErrorOk;
	}

	// ============================================================================
	// [mathpresso::Context - Construction / Destruction]
	// ============================================================================

	Context::Context()
		: _symbols(std::make_shared<Symbols>()),
		_parent(),
		_children({}),
		_isGlobal(true)
	{
	}

	Context::Context(const Context& other)
		:_symbols(other._symbols),
		_parent(other._parent),
		_children(other._children),
		_isGlobal(other._isGlobal)
	{
	}

	// ============================================================================
	// [mathpresso::Context - Interface]
	// ============================================================================

	Error Context::addBuiltIns(void)
	{
		addBuiltinMpObjects(this);
		return ErrorCode::kErrorOk;
	}

	Error Context::addConstant(const std::string &name, double value)
	{
		auto shared_sym(std::make_shared<AstSymbol>(name, AstSymbolType::kAstSymbolVariable, AstScopeType::kAstScopeGlobal));

		shared_sym->setValue(value);
		shared_sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsDeclared | AstSymbolFlags::kAstSymbolIsReadOnly | AstSymbolFlags::kAstSymbolIsAssigned);

		_symbols->add(name, shared_sym);

		return ErrorCode::kErrorOk;
	}

	Error Context::addConstant(const std::string &name, std::complex<double> value)
	{
		auto shared_sym(std::make_shared<AstSymbol>(name, AstSymbolType::kAstSymbolVariable, AstScopeType::kAstScopeGlobal));

		shared_sym->setValue(value);
		shared_sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsDeclared | AstSymbolFlags::kAstSymbolIsReadOnly | AstSymbolFlags::kAstSymbolIsAssigned | AstSymbolFlags::kAstSymbolIsComplex);

		_symbols->add(name, shared_sym);

		return ErrorCode::kErrorOk;
	}

	Error Context::addVariable(const std::string &name, int offset, unsigned int flags)
	{
		auto shared_sym(std::make_shared<AstSymbol>(name, AstSymbolType::kAstSymbolVariable, AstScopeType::kAstScopeGlobal));

		shared_sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsDeclared);
		if (flags & VariableFlags::kVariableCplx)
		{
			shared_sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsComplex);
		}

		shared_sym->setVarSlotId(InternalConsts::kInvalidSlot);
		shared_sym->setVarOffset(offset);

		if (flags & VariableFlags::kVariableRO)
		{
			shared_sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsReadOnly);
		}

		_symbols->add(name, shared_sym);
		return ErrorCode::kErrorOk;
	}

	Error Context::addObject(const std::string &name, std::shared_ptr<MpOperation> obj)
	{
		_symbols->add(name, obj);
		return ErrorCode::kErrorOk;
	}

	Error Context::listSymbols(std::vector<std::string> &syms)
	{
		syms = _symbols->names();
		return ErrorCode::kErrorOk;
	}

	std::vector<std::shared_ptr<AstSymbol>> Context::getVariables()
	{
		std::vector<std::shared_ptr<AstSymbol>> variables, tmp;
		std::shared_ptr<Context> ctx = shared_from_this();
		do
		{
			tmp = ctx->_symbols->getVariables();
			variables.insert(variables.end(), tmp.begin(), tmp.end());
		} while (ctx = ctx->getParent());

		return variables;
	}

	inline std::shared_ptr<Context> Context::getChild(const std::string & name)
	{
		try
		{
			return _children.at(name);
		}
		catch (std::out_of_range)
		{
			return nullptr;
		}
	}

	Error Context::setParent(std::shared_ptr<Context> ctx)
	{
		_parent = ctx;
		return ErrorCode::kErrorOk;
	}

	Error Context::addChild(const std::string & name, std::shared_ptr<Context> ctx)
	{
		auto ret = _children.emplace(name, ctx);
		if (!ret.second)
			return ErrorCode::kErrorInvalidArgument;

		return ctx->setParent(shared_from_this());
	}

	Error Context::delSymbol(const std::string &name)
	{
		_symbols->remove(name);
		return ErrorCode::kErrorOk;
	}

	// ============================================================================
	// [mathpresso::Expression - Construction / Destruction]
	// ============================================================================

	Expression::Expression() : _func(mpDummyFunc)
	{
	}

	Expression::~Expression()
	{
		reset();
	}

	// ============================================================================
	// [mathpresso::Expression - Interface]
	// ============================================================================

	Error Expression::compile(std::shared_ptr<Context> ctx, const std::string & body, unsigned int options, OutputLog* log)
	{
		// Init options first.
		options &= _kOptionsMask;

		if (log != nullptr)
			options |= InternalOptions::kInternalOptionLog;
		else
			options &= ~(kOptionVerbose | kOptionDebugAst | kOptionDebugAsm);

		StringBuilderTmp<512> sbTmp;

		// Initialize AST.

		// this AstBuilder will hold the AST and the symbols defined by assignment within the expression.
		// Every other (global) Variable will be hold within ctx._d!
		std::shared_ptr<AstBuilder> ast(std::make_shared<AstBuilder>());
		MATHPRESSO_PROPAGATE(ast->initProgramScope());

		// Setup basic data structures used during parsing and compilation.
		ErrorReporter errorReporter(body, options, log);

		if (options & Options::kOptionVerbose)
		{
			log->log(OutputLog::kMessageWarning, 0, 0, body);
		}

		// create shadowContext and add ctx as parent. here all expression-local symbols will be stored.
		std::shared_ptr<Context> shadowContext(std::make_shared<Context>());
		shadowContext->markShadow();
		shadowContext->setParent(ctx);

		// Parse the expression into AST.
		{
			MATHPRESSO_PROPAGATE(Parser(ast, &errorReporter, body, shadowContext).parseProgram(ast->getProgramNode()));
		}

		if (options & kOptionDebugAst)
		{
			ast->dump(sbTmp, ctx->_symbols);
			log->log(OutputLog::kMessageAstInitial, 0, 0, sbTmp.getData(), sbTmp.getLength());
			sbTmp.clear();
		}

		// Perform basic optimizations at AST level.
		{
			MATHPRESSO_PROPAGATE(AstOptimizer(ast, &errorReporter, shadowContext).onProgram(ast->getProgramNode()));
		}

		if (options & kOptionDebugAst)
		{
			ast->dump(sbTmp, ctx->_symbols);
			log->log(OutputLog::kMessageAstFinal, 0, 0, sbTmp.getData(), sbTmp.getLength());
			sbTmp.clear();
		}

		_isComplex = ast->_programNode->returnsComplex();

		// Compile the function to machine code.
		reset();
		CompiledFunc fn = mpCompileFunction(ast, options, log, ctx, _isComplex);

		if (fn == nullptr)
			return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorNoMemory);
		_func = fn;

		return ErrorCode::kErrorOk;
	}

	bool Expression::isCompiled() const
	{
		return _func != mpDummyFunc;
	}

	void Expression::reset()
	{
		// Allocated by a JIT memory manager, free it.
		if (_func != mpDummyFunc)
		{
			mpFreeFunction((void*)_func);
			_func = mpDummyFunc;
		}
	}

	// ============================================================================
	// [mathpresso::OutputLog - Construction / Destruction]
	// ============================================================================

	OutputLog::OutputLog()
	{
	}
	OutputLog::~OutputLog()
	{
	}

	// ============================================================================
	// [mathpresso::ErrorReporter - Interface]
	// ============================================================================

	void ErrorReporter::getLineAndColumn(uint32_t position, uint32_t& line, uint32_t& column)
	{
		// Shouldn't happen, but be defensive.
		if (static_cast<size_t>(position) >= _len)
		{
			line = 0;
			column = 0;
			return;
		}

		const char* pStart = _body;
		const char* p = pStart + position;

		uint32_t x = 0;
		uint32_t y = 1;

		// Find column.
		while (p[0] != '\n')
		{
			x++;
			if (p == pStart)
				break;
			p--;
		}

		// Find line.
		while (p != pStart)
		{
			y += p[0] == '\n';
			p--;
		}

		line = y;
		column = x;
	}

	void ErrorReporter::onWarning(uint32_t position, const char* fmt, ...)
	{
		if (reportsWarnings())
		{
			StringBuilderTmp<256> sb;

			va_list ap;
			va_start(ap, fmt);

			sb.appendFormatVA(fmt, ap);

			va_end(ap);
			onWarning(position, sb);
		}
	}

	void ErrorReporter::onWarning(uint32_t position, const StringBuilder& msg)
	{
		if (reportsWarnings())
		{
			uint32_t line, column;
			getLineAndColumn(position, line, column);
			_log->log(OutputLog::kMessageWarning, line, column, msg.getData(), msg.getLength());
		}
	}

	void ErrorReporter::onWarning(uint32_t position, const std::string & msg)
	{
		if (reportsWarnings())
		{
			uint32_t line, column;
			getLineAndColumn(position, line, column);
			_log->log(OutputLog::kMessageWarning, line, column, msg.c_str(), msg.length());
		}
	}

	Error ErrorReporter::onError(Error error, uint32_t position, const char* fmt, ...)
	{
		if (reportsErrors())
		{
			StringBuilderTmp<256> sb;

			va_list ap;
			va_start(ap, fmt);

			sb.appendFormatVA(fmt, ap);

			va_end(ap);
			return onError(error, position, sb);
		}
		else
		{
			return MATHPRESSO_TRACE_ERROR(error);
		}
	}

	Error ErrorReporter::onError(Error error, uint32_t position, const StringBuilder& msg)
	{
		if (reportsErrors())
		{
			uint32_t line, column;
			getLineAndColumn(position, line, column);
			_log->log(OutputLog::kMessageError, line, column, msg.getData(), msg.getLength());
		}

		return MATHPRESSO_TRACE_ERROR(error);
	}

	Error ErrorReporter::onError(Error error, uint32_t position, const std::string & msg)
	{
		if (reportsErrors())
		{
			uint32_t line, column;
			getLineAndColumn(position, line, column);
			_log->log(OutputLog::kMessageError, line, column, msg.c_str(), msg.length());
		}

		return MATHPRESSO_TRACE_ERROR(error);
	}


	std::string Symbols::name(const std::shared_ptr<MpOperation> ptr) const
	{
		for (auto &pn : _operations)
		{
			for (auto &p : pn.second)
			{
				if (p == ptr)
					return pn.first;
			}
		}
		return "<unknown>";
	}

	std::shared_ptr<MpOperation> Symbols::findFunction(const std::string &name, size_t numArgs) const
	{
		auto ps = findFunction(name);
		for (auto &p : ps)
		{
			if (p->nargs() == numArgs)
				return p;
		}

		return nullptr;
	}

	std::vector<Symbols::op_ptr_type> Symbols::findFunction(const std::string &name) const
	{
		auto it = _operations.find(name);
		if (it == _operations.end())
			return{};
		else
			return it->second;
	}

	Symbols::op_ptr_type Symbols::findFunction(const std::string & name, size_t nargs, bool paramsAreComplex) const
	{
		auto it = _operations.find(name);

		if (it == _operations.end())
			return nullptr;

		Symbols::op_ptr_type weakFit = nullptr;

		// use a reverse iterator, as overrides are added after the 'original'.
		for (auto p = it->second.rbegin(); p != it->second.rend(); p++)
		{
			// Need the right number of arguments
			if (p->get()->nargs() == nargs)
			{
				if (paramsAreComplex)
				{
					if (p->get()->signature().areParams(Signature::type::complex))
						return *p;
				}
				else
				{
					if (p->get()->signature().areParams(Signature::type::real))
						return *p;
					else if (!weakFit)
					{
						weakFit = *p;
					}
				}
			}
		}
		return weakFit;
	}

	Symbols::var_ptr_type Symbols::findVariable(const std::string & name) const
	{
		try
		{
			return _variables.at(name);
		}
		catch (std::out_of_range)
		{
			return nullptr;
		}
	}

	void Symbols::add(const std::string &name, Symbols::op_ptr_type obj)
	{
		auto syms = _operations[name];
		for (auto p : syms)
		{
			if (p->nargs() == obj->nargs())
			{
				if (p->precedence() != obj->precedence()
					|| (p->flags() & MpOperation::Flags::RighttoLeft) != (obj->flags() & MpOperation::Flags::RighttoLeft))
				{
					throw std::runtime_error("unable to add function.");
				}
			}
		}


		_operations[name].push_back(obj);
	}

	void Symbols::add(const std::string & name, var_ptr_type obj)
	{
		if (_variables.find(name) == _variables.end())
		{
			_variables[name] = obj;
		}
	}

	void Symbols::remove(const std::string &name)
	{
		auto it = _operations.find(name);
		if (it != _operations.end())
			_operations.erase(it);
	}

	std::vector<std::string> Symbols::names() const
	{
		std::vector<std::string> names;
		for (auto &pn : _operations)
		{
			for (auto &p : pn.second)
				names.push_back(pn.first + " (" + p->signature().to_string() + ")");
		}

		for (auto p : _variables)
		{
			names.push_back(p.first + (p.second->hasSymbolFlag(AstSymbolFlags::kAstSymbolIsComplex) ? " (complex)" : " (real)"));
		}

		return names;
	}

	void Symbols::clear()
	{
		_operations.clear();
		_variables.clear();
	}

	std::vector<std::shared_ptr<AstSymbol>> Symbols::getVariables()
	{
		std::vector<std::shared_ptr<AstSymbol>> out;
		for (auto p : _variables)
		{
			out.push_back(p.second);
		}
		return out;
	}

	namespace resolver
	{
		std::shared_ptr<MpOperation> resolveFunction(ContextPtr ctx, const std::string & name, size_t numargs, bool takesComplex)
		{
			std::shared_ptr<MpOperation> function;
			auto splitName(separateName(name));
			do
			{
				ContextPtr tmpCtx(ctx);
				for (size_t i = 0; tmpCtx && i < splitName.size() - 1; i++)
				{
					tmpCtx = tmpCtx->getChild(splitName[i]);
				}

				if (tmpCtx)
				{
					function = tmpCtx->_symbols->findFunction(splitName.back(), numargs, takesComplex);
				}
			} while (!function && (ctx = ctx->getParent()) != nullptr);


			return function;
		}

		std::vector<std::shared_ptr<MpOperation>> resolveFunction(ContextPtr ctx, const std::string & name)
		{
			std::vector<std::shared_ptr<MpOperation>> functions, tmp;
			auto splitName(separateName(name));
			do
			{
				ContextPtr tmpCtx(ctx);
				for (size_t i = 0; tmpCtx && i < splitName.size() - 1; i++)
				{
					tmpCtx = tmpCtx->getChild(splitName[i]);
				}

				if (tmpCtx)
				{
					tmp = tmpCtx->_symbols->findFunction(splitName.back());
					functions.insert(functions.begin(), tmp.begin(), tmp.end());
				}
			} while (ctx = ctx->getParent());

			return functions;
		}

		std::shared_ptr<MpOperation> resolveFunction(ContextPtr ctx, const std::string & name, size_t numargs)
		{
			std::shared_ptr<MpOperation> function;
			auto splitName(separateName(name));
			do
			{
				ContextPtr tmpCtx(ctx);
				for (size_t i = 0; tmpCtx && i < splitName.size() - 1; i++)
				{
					tmpCtx = tmpCtx->getChild(splitName[i]);
				}

				if (tmpCtx)
				{
					function = tmpCtx->_symbols->findFunction(splitName.back(), numargs);
				}
			} while (!function && (ctx = ctx->getParent()) != nullptr);

			return function;
		}

		std::shared_ptr<AstSymbol> resolveVariable(ContextPtr ctx, const std::string & name, ContextPtr * contextOut)
		{
			std::shared_ptr<AstSymbol> symbol;

			auto splitName(separateName(name));

			do
			{
				ContextPtr tmpCtx(ctx);
				for (size_t i = 0; tmpCtx && i < splitName.size() - 1; i++)
				{
					tmpCtx = tmpCtx->getChild(splitName[i]);
				}

				if (tmpCtx)
				{
					symbol = tmpCtx->_symbols->findVariable(splitName.back());
				}

			} while (symbol == nullptr && (ctx = ctx->getParent()) != nullptr);

			if (contextOut != nullptr)
				*contextOut = ctx;

			return symbol;
		}

		std::vector<std::string> separateName(std::string name)
		{
			std::vector<std::string> ret;

			std::istringstream iss(name);
			std::string token;
			while (std::getline(iss, token, '.'))
			{
				ret.push_back(token);
			}


			return ret;
		}
	}

} // mathpresso namespace
