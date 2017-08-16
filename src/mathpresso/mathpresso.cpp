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

	MATHPRESSO_NOAPI Error mpTraceError(Error error)
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
	// [mathpresso::ContextInternalImpl]
	// ============================================================================

	//! \internal
	//!
	//! Null context data.
	static const ContextImpl mpContextNull = { 0 };

	//! \internal
	//!
	//! Internal context data.
	struct ContextInternalImpl : public ContextImpl
	{
		MATHPRESSO_INLINE ContextInternalImpl()
			: _zone(32768 - Zone::kZoneOverhead),
			_heap(&_zone),
			_builder(&_heap),
			_scope(&_builder, static_cast<AstScope*>(NULL), AstScopeType::kAstScopeGlobal)
		{
			mpAtomicSet(&_refCount, 1);
		}
		MATHPRESSO_INLINE ~ContextInternalImpl() {}

		Zone _zone;
		ZoneHeap _heap;
		AstBuilder _builder;
		AstScope _scope;
	};

	static MATHPRESSO_INLINE ContextImpl* mpContextAddRef(ContextImpl* d)
	{
		if (d != &mpContextNull)
			mpAtomicInc(&d->_refCount);
		return d;
	}

	static MATHPRESSO_INLINE void mpContextRelease(ContextImpl* d)
	{
		if (d != &mpContextNull && !mpAtomicDec(&d->_refCount))
			delete static_cast<ContextInternalImpl*>(d);
	}

	static ContextImpl* mpContextClone(ContextImpl* otherD_)
	{
		ContextInternalImpl* d = new(std::nothrow) ContextInternalImpl();
		if (MATHPRESSO_UNLIKELY(d == nullptr)) return nullptr;

		if (otherD_ != &mpContextNull)
		{
			ContextInternalImpl* otherD = static_cast<ContextInternalImpl*>(otherD_);
			AstSymbolHashIterator it(otherD->_scope._symbols);

			while (it.has())
			{
				AstSymbol* sym = it.get();

				std::string name(sym->getName(), sym->getLength());
				uint32_t hVal = sym->getHVal();
				uint32_t type = sym->getSymbolType();

				AstSymbol* clonedSym = d->_builder.newSymbol(name, hVal, type, otherD->_scope.getScopeType());
				if (MATHPRESSO_UNLIKELY(clonedSym == nullptr))
				{
					delete d;
					return nullptr;
				}

				clonedSym->setSymbolFlag(sym->getSymbolFlags());
				switch (type)
				{
					case AstSymbolType::kAstSymbolVariable:
						clonedSym->setVarSlotId(sym->getVarSlotId());
						clonedSym->setVarOffset(sym->getVarOffset());
						clonedSym->setValue(sym->getValueComp());
						break;

					case AstSymbolType::kAstSymbolIntrinsic:
					case AstSymbolType::kAstSymbolFunction:
						break;

					default:
						MATHPRESSO_ASSERT_NOT_REACHED();
				}

				d->_scope.putSymbol(clonedSym);
				it.next();
			}
		}

		return d;
	}

	static Error mpContextMutable(Context* self, ContextInternalImpl** out)
	{
		ContextImpl* d = self->_d;

		if (d != &mpContextNull && mpAtomicGet(&d->_refCount) == 1)
		{
			*out = static_cast<ContextInternalImpl*>(d);
			return ErrorCode::kErrorOk;
		}
		else
		{
			d = mpContextClone(d);
			if (MATHPRESSO_UNLIKELY(d == nullptr))
				return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorNoMemory);

			mpContextRelease(
				mpAtomicSetXchgT<ContextImpl*>(
					&self->_d, d));

			*out = static_cast<ContextInternalImpl*>(d);
			return ErrorCode::kErrorOk;
		}
	}


	// ============================================================================
	// [mathpresso::Context - Copy / Reset]
	// ============================================================================

	Error Context::reset()
	{
		mpContextRelease(
			mpAtomicSetXchgT<ContextImpl*>(
				&_d, const_cast<ContextImpl*>(&mpContextNull)));

		return ErrorCode::kErrorOk;
	}

	Context& Context::operator=(const Context& other)
	{
		mpContextRelease(
			mpAtomicSetXchgT<ContextImpl*>(
				&_d, mpContextAddRef(other._d)));
		return *this;
	}


	// ============================================================================
	// [mathpresso::Context - Construction / Destruction]
	// ============================================================================

	Context::Context()
		: _d(const_cast<ContextImpl*>(&mpContextNull)),
		_ops()
	{
	}

	Context::Context(const Context& other)
		: _d(mpContextAddRef(other._d)),
		_ops(other._ops)
	{
	}

	Context::~Context()
	{
		mpContextRelease(_d);
	}


	// ============================================================================
	// [mathpresso::Context - Interface]
	// ============================================================================

	Error Context::addBuiltIns(void)
	{
		addBuiltinMpObjects(this);
		return ErrorCode::kErrorOk;
	}

	Error Context::addSymbol(AstSymbol* &sym, const char * name, int type)
	{
		ContextInternalImpl* d;
		MATHPRESSO_PROPAGATE(mpContextMutable(this, &d));

		size_t nlen = strlen(name);
		uint32_t hVal = HashUtils::hashString(name, nlen);
		sym = d->_scope.getSymbol(std::string(name, nlen), hVal);
		if (sym != nullptr)
			return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorSymbolAlreadyExists);

		sym = d->_builder.newSymbol(std::string(name, nlen), hVal, type, AstScopeType::kAstScopeGlobal);
		if (sym == nullptr)
			return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorNoMemory);
		d->_scope.putSymbol(sym);

		return ErrorCode::kErrorOk;
	}

	Error Context::addConstant(const char* name, double value)
	{
		AstSymbol* sym;
		MATHPRESSO_PROPAGATE(addSymbol(sym, name, AstSymbolType::kAstSymbolVariable));

		sym->setValue(value);
		sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsDeclared | AstSymbolFlags::kAstSymbolIsReadOnly | AstSymbolFlags::kAstSymbolIsAssigned);

		return ErrorCode::kErrorOk;
	}

	Error Context::addConstant(const char* name, std::complex<double> value)
	{
		AstSymbol* sym;
		MATHPRESSO_PROPAGATE(addSymbol(sym, name, AstSymbolType::kAstSymbolVariable));

		sym->setValue(value);
		sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsDeclared | AstSymbolFlags::kAstSymbolIsReadOnly | AstSymbolFlags::kAstSymbolIsAssigned | AstSymbolFlags::kAstSymbolIsComplex);

		return ErrorCode::kErrorOk;
	}

	Error Context::addVariable(const char* name, int offset, unsigned int flags)
	{
		AstSymbol* sym;
		MATHPRESSO_PROPAGATE(addSymbol(sym, name, AstSymbolType::kAstSymbolVariable));

		sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsDeclared);
		if (flags & VariableFlags::kVariableCplx)
			sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsComplex);

		sym->setVarSlotId(InternalConsts::kInvalidSlot);
		sym->setVarOffset(offset);

		if (flags & VariableFlags::kVariableRO)
			sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsReadOnly);

		return ErrorCode::kErrorOk;
	}

	Error Context::addObject(std::string name, std::shared_ptr<MpOperation> obj)
	{
		AstSymbol * sym;
		Error e = addSymbol(sym, name.c_str(), AstSymbolType::kAstSymbolFunction);
		if (e != ErrorCode::kErrorOk)
		{
			if (e != ErrorCode::kErrorSymbolAlreadyExists || sym->getSymbolType() != AstSymbolType::kAstSymbolFunction)
				return e;
		}
		else
		{
			sym->setSymbolFlag(kAstSymbolIsDeclared);
		}

		_ops.add(name, obj);

		return ErrorCode::kErrorOk;
	}

	Error Context::listSymbols(std::vector<std::string> &syms)
	{
		syms = _ops.names();

		ContextInternalImpl* d;
		MATHPRESSO_PROPAGATE(mpContextMutable(this, &d));

		HashIterator<std::string, AstSymbol> it(d->_scope.getSymbols());
		do
		{
			if (it.get()->getSymbolType() == AstSymbolType::kAstSymbolVariable)
			{
				syms.push_back(it.get()->getName());
			}
		} while (it.next());

		return ErrorCode::kErrorOk;
	}

	Error Context::delSymbol(const char* name)
	{
		ContextInternalImpl* d;
		MATHPRESSO_PROPAGATE(mpContextMutable(this, &d));

		size_t nlen = strlen(name);
		uint32_t hVal = HashUtils::hashString(name, nlen);

		AstSymbol* sym = d->_scope.getSymbol(std::string(name, nlen), hVal);
		if (sym == nullptr)
			return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorSymbolNotFound);

		_ops.remove(name, sym->getLength());
		d->_builder.deleteSymbol(d->_scope.removeSymbol(sym));

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

	Error Expression::compile(const Context& ctx, const char* body, unsigned int options, OutputLog* log)
	{
		// Init options first.
		options &= _kOptionsMask;

		if (log != nullptr)
			options |= InternalOptions::kInternalOptionLog;
		else
			options &= ~(kOptionVerbose | kOptionDebugAst | kOptionDebugAsm);

		Zone zone(32768 - Zone::kZoneOverhead);
		ZoneHeap heap(&zone);
		StringBuilderTmp<512> sbTmp;

		// Initialize AST.
		AstBuilder ast(&heap);
		MATHPRESSO_PROPAGATE(ast.initProgramScope());

		ContextImpl* d = ctx._d;
		if (d != &mpContextNull)
			ast.getRootScope()->shadowContextScope(&static_cast<ContextInternalImpl*>(d)->_scope);

		// Setup basic data structures used during parsing and compilation.
		size_t len = ::strlen(body);
		ErrorReporter errorReporter(body, len, options, log);

		// Parse the expression into AST.
		{
			MATHPRESSO_PROPAGATE(Parser(&ast, &errorReporter, body, len, &ctx._ops).parseProgram(ast.getProgramNode())); 
		}

		if (options & kOptionDebugAst)
		{
			ast.dump(sbTmp, &ctx._ops);
			log->log(OutputLog::kMessageAstInitial, 0, 0, sbTmp.getData(), sbTmp.getLength());
			sbTmp.clear();
		}

		// Perform basic optimizations at AST level.
		{
			MATHPRESSO_PROPAGATE(AstOptimizer(&ast, &errorReporter, &ctx._ops).onProgram(ast.getProgramNode()));
		}

		if (options & kOptionDebugAst)
		{
			ast.dump(sbTmp, &ctx._ops);
			log->log(OutputLog::kMessageAstFinal, 0, 0, sbTmp.getData(), sbTmp.getLength());
			sbTmp.clear();
		}

		_isComplex = ast._programNode->returnsComplex();

		// Compile the function to machine code.
		reset();
		CompiledFunc fn = mpCompileFunction(&ast, options, log, &ctx._ops, _isComplex);

		if (fn == nullptr)
			return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorNoMemory);
		_func = fn;

		//ast.getProgramNode();
		//ast.deleteNode(ast.getProgramNode());

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

	std::string Operations::name(const MpOperation * ptr) const
	{
		auto test = std::find_if(_symbols.begin(), _symbols.end(), [&](const std::pair<std::pair<std::string, size_t>, std::shared_ptr<MpOperation>> &pair)
		{
			return pair.second.get() == ptr;
		});

		if (test != _symbols.end())
			return test->first.first;
		else
			return "<unknown>";
	}

	std::string Operations::name(const std::shared_ptr<MpOperation> ptr) const
	{
		return name(ptr.get());
	}


	std::shared_ptr<MpOperation> Operations::find(const std::string &name, size_t numArgs) const
	{
		auto ret = _symbols.find(std::make_pair(name, numArgs));
		if (ret == _symbols.end())
		{
			 return nullptr;
		}
		else
		{
			return ret->second;
		}
	}

	std::vector<Operations::op_ptr_type> Operations::find(const std::string &name) const
	{
		std::vector<op_ptr_type> ret;

		for (auto p : _symbols)
		{
			if (p.first.first == name)
				ret.push_back(p.second);
		}

		return ret;
	}

	Operations::op_ptr_type Operations::find(const std::string & name, size_t nargs, bool paramsAreComplex) const
	{
		Operations::op_ptr_type weakFit = nullptr;
		for (auto p : _symbols)
		{
			if (p.first.first == name && p.first.second == nargs)
			{
				if (paramsAreComplex)
				{
					if (p.second->signature().areParams(Signature::type::complex))
						return p.second;
				}
				else
				{
					if (p.second->signature().areParams(Signature::type::real))
						return p.second;
					else if (!weakFit)
					{
						weakFit = p.second;
					}
				}
			}
		}
		return weakFit;
	}

	void Operations::add(const std::string &name, std::shared_ptr<MpOperation> obj)
	{
		_symbols.emplace(std::make_pair(name, obj->nargs()), obj);
	}

	void Operations::remove(const std::string &name, size_t numArgs)
	{
		auto it = _symbols.find(std::make_pair(name, numArgs));
		if (it != _symbols.end())
			_symbols.erase(it);
	}

	std::vector<std::string> Operations::names() const
	{
		std::vector<std::string> names;
		for (auto p : _symbols)
		{
			names.push_back(p.first.first + " (" + p.second->signature().to_string() + ")");
		}
		return names;
	}

} // mathpresso namespace
