// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include "./mathpresso_p.h"
#include "./mpast_p.h"
#include "./mpatomic_p.h"
#include "./mpcompiler_p.h"
#include "./mpeval_p.h"
#include "./mpoptimizer_p.h"
#include "./mpparser_p.h"
#include "./mptokenizer_p.h"

#include <math.h>
#include <string.h>

namespace mathpresso {

// ============================================================================
// [mathpresso::mpOpInfo]
// ============================================================================

// Operator information, precedence and association. The table is mostly based
// on the C-language standard, but also adjusted to support MATHPRESSO specific
// operators and rules. However, the associativity and precedence should be
// fully compatible with C.

#define LTR 0
#define RTL kOpFlagRightToLeft
#define CtoC (kOpFlagComplexToComplex)
#define DtoD (kOpFlagRealToReal)
#define CandD (DtoD | CtoC)
#define CtoD (kOpFlagComplexToReal)
#define DtoC (kOpFlagRealToComplex)
	std::vector<std::pair<std::string, OpInfo>> _symbols = {
		{ "_none_", OpInfo("_none_", kOpNone, 0, LTR) },
		{ "-", OpInfo("-", kOpNeg, 3, RTL | CandD | kOpFlagArithmetic | kOpFlagUnary) },
		{ "!", OpInfo("!", kOpNot, 3, RTL | DtoD | kOpFlagCondition | kOpFlagUnary) },
		{ "=", OpInfo("=", kOpAssign, 15, RTL | CandD | kOpFlagAssign | kOpFlagBinary) },
		{ "==", OpInfo("==", kOpEq, 9, LTR | CandD | kOpFlagCondition | kOpFlagBinary) },
		{ "!=", OpInfo("!=", kOpNe, 9, LTR | CandD | kOpFlagCondition | kOpFlagBinary) },
		{ "<", OpInfo("<", kOpLt, 8, LTR | DtoD | kOpFlagCondition | kOpFlagBinary) },
		{ "<=", OpInfo("<=", kOpLe, 8, LTR | DtoD | kOpFlagCondition | kOpFlagBinary) },
		{ ">", OpInfo(">", kOpGt, 8, LTR | DtoD | kOpFlagCondition | kOpFlagBinary) },
		{ "<=", OpInfo("<=", kOpGe, 8, LTR | DtoD | kOpFlagCondition | kOpFlagBinary) },
		{ "+", OpInfo("+", kOpAdd, 6, LTR | CandD | kOpFlagArithmetic | kOpFlagNopIfZero | kOpFlagBinary, mpAddOptD, mpAddOptC, nullptr, nullptr, compileAddD, compileAddC) },
		{ "-", OpInfo("-", kOpSub, 6, LTR | CandD | kOpFlagArithmetic | kOpFlagNopIfRZero | kOpFlagBinary) },
		{ "*", OpInfo("*", kOpMul, 5, LTR | CandD | kOpFlagArithmetic | kOpFlagNopIfOne | kOpFlagBinary) },
		{ "/", OpInfo("/", kOpDiv, 5, LTR | CandD | kOpFlagArithmetic | kOpFlagNopIfROne | kOpFlagBinary) },
		{ "?", OpInfo("?", kOpQMark, 15, RTL | kOpFlagTernary) },
		{ ":", OpInfo(":", kOpColon, 15, RTL | kOpFlagTernary) },
		{ "%", OpInfo("%", kOpMod, 5, LTR | DtoD | kOpFlagBinary | kOpFlagIntrinsic, mpMod, nullptr) },
		{ "isnan", OpInfo("isnan", kOpIsNan, 0, LTR | DtoD | kOpFlagCondition | kOpFlagUnary | kOpFlagIntrinsic, mpIsNan, nullptr) },
		{ "isinf", OpInfo("isinf", kOpIsInf, 0, LTR | DtoD | kOpFlagCondition | kOpFlagUnary | kOpFlagIntrinsic, mpIsInf, nullptr) },
		{ "isfinite", OpInfo("isfinite", kOpIsFinite, 0, LTR | DtoD | kOpFlagCondition | kOpFlagUnary | kOpFlagIntrinsic, mpIsFinite, nullptr) },
		{ "signbit", OpInfo("signbit", kOpSignBit, 0, LTR | DtoD | kOpFlagCondition | kOpFlagUnary | kOpFlagIntrinsic, mpSignBit, nullptr) },
		{ "round", OpInfo("round", kOpRound, 0, LTR | DtoD | kOpFlagRounding | kOpFlagUnary | kOpFlagIntrinsic, mpRound, nullptr) },
		{ "roundeven", OpInfo("roundeven", kOpRoundEven, 0, LTR | DtoD | kOpFlagRounding | kOpFlagUnary | kOpFlagIntrinsic, mpRoundEven, nullptr) },
		{ "trunc", OpInfo("trunc", kOpTrunc, 0, LTR | DtoD | kOpFlagRounding | kOpFlagUnary | kOpFlagIntrinsic, mpTrunc, nullptr) },
		{ "floor", OpInfo("floor", kOpFloor, 0, LTR | DtoD | kOpFlagRounding | kOpFlagUnary | kOpFlagIntrinsic, mpFloor, nullptr) },
		{ "ceil", OpInfo("ceil", kOpCeil, 0, LTR | DtoD | kOpFlagRounding | kOpFlagUnary | kOpFlagIntrinsic, mpCeil, nullptr) },
		{ "abs", OpInfo("abs", kOpAbs, 0, LTR | DtoD | CtoD | kOpFlagUnary | kOpFlagIntrinsic, mpAbs, nullptr, mpAbsC, nullptr) },
		{ "exp", OpInfo("exp", kOpExp, 0, LTR | CandD | kOpFlagUnary | kOpFlagIntrinsic, mpExp, mpFuncCtoC1<std::exp>) },
		{ "log", OpInfo("log", kOpLog, 0, LTR | CandD | kOpFlagUnary | kOpFlagIntrinsic, mpLog, mpFuncCtoC1<std::log>) },
		{ "log2", OpInfo("log2", kOpLog2, 0, LTR | CandD | kOpFlagUnary | kOpFlagIntrinsic, mpLog2, mpLog2C) },
		{ "log10", OpInfo("log10", kOpLog10, 0, LTR | CandD | kOpFlagUnary | kOpFlagIntrinsic, mpLog10, mpFuncCtoC1<std::log10>) },
		{ "sqrt", OpInfo("sqrt", kOpSqrt, 0, LTR | DtoD | kOpFlagUnary | kOpFlagIntrinsic, mpSqrt, nullptr) },
		{ "frac", OpInfo("frac", kOpFrac, 0, LTR | DtoD | kOpFlagUnary | kOpFlagIntrinsic, mpFrac, nullptr) },
		{ "recip", OpInfo("recip", kOpRecip, 0, LTR | CandD | kOpFlagUnary | kOpFlagIntrinsic, mpRecip, mpRecipC) },
		{ "sin", OpInfo("sin", kOpSin, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpSin, mpFuncCtoC1<std::sin>) },
		{ "cos", OpInfo("cos", kOpCos, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpCos, mpFuncCtoC1<std::cos>) },
		{ "tan", OpInfo("tan", kOpTan, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpTan, mpFuncCtoC1<std::tan>) },
		{ "sinh", OpInfo("sinh", kOpSinh, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpSinh, mpFuncCtoC1<std::sinh>) },
		{ "cosh", OpInfo("cosh", kOpCosh, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpCosh, mpFuncCtoC1<std::cosh>) },
		{ "tanh", OpInfo("tanh", kOpTanh, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpTanh, mpFuncCtoC1<std::tanh>) },
		{ "asin", OpInfo("asin", kOpAsin, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpAsin, mpFuncCtoC1<std::asin>) },
		{ "acos", OpInfo("acos", kOpAcos, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpAcos, mpFuncCtoC1<std::acos>) },
		{ "atan", OpInfo("atan", kOpAtan, 0, LTR | CandD | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagIntrinsic, mpAtan, mpFuncCtoC1<std::atan>) },
		{ "avg", OpInfo("avg", kOpAvg, 0, LTR | CandD | kOpFlagBinary | kOpFlagIntrinsic, mpAvg, mpAvgC) },
		{ "min", OpInfo("min", kOpMin, 0, LTR | DtoD | kOpFlagBinary | kOpFlagIntrinsic, mpMin<double>, nullptr) },
		{ "max", OpInfo("max", kOpMax, 0, LTR | DtoD | kOpFlagBinary | kOpFlagIntrinsic, mpMax<double>, nullptr) },
		{ "pow", OpInfo("pow", kOpPow, 0, LTR | CandD | kOpFlagBinary | kOpFlagNopIfROne | kOpFlagIntrinsic, mpPow, mpFuncCtoC2<std::pow>) },
		{ "atan2", OpInfo("atan2", kOpAtan2, 0, LTR | DtoD | kOpFlagTrigonometric | kOpFlagBinary | kOpFlagIntrinsic, mpAtan2, nullptr) },
		{ "hypot", OpInfo("hypot", kOpHypot, 0, LTR | DtoD | kOpFlagTrigonometric | kOpFlagBinary | kOpFlagIntrinsic, mpHypot, nullptr) },
		{ "copysign", OpInfo("copysign", kOpCopySign, 0, LTR | DtoD | kOpFlagBinary | kOpFlagIntrinsic, mpCopySign, nullptr) },
		{ "real", OpInfo("real", kOpReal, 0, LTR | CtoD | kOpFlagUnary | kOpFlagIntrinsic, nullptr, nullptr, mpGetReal, nullptr) },
		{ "imag", OpInfo("imag", kOpImag, 0, LTR | CtoD | kOpFlagUnary | kOpFlagIntrinsic, nullptr, nullptr, mpGetImag, nullptr) },
		{ "conjug", OpInfo("conjug", kOpConjug, 0, LTR | CtoC | kOpFlagUnary | kOpFlagIntrinsic, nullptr, mpFuncCtoC1<std::conj>) },
		{ "sqrtC", OpInfo("sqrtC", kOpSqrtC, 0, LTR | CtoC | kOpFlagUnary | kOpFlagIntrinsic, nullptr, mpFuncCtoC1<std::sqrt>) }
	};
#undef DtoC
#undef CtoD
#undef CandD
#undef DtoD
#undef CtoC
#undef RTL
#undef LTR

// ============================================================================
// [mathpresso::mpAssertionFailure]
// ============================================================================

void mpAssertionFailure(const char* file, int line, const char* msg) {
  fprintf(stderr,
    "[mathpresso] Assertion failed at %s (line %d):\n"
    "[mathpresso] %s\n", file, line, msg);

  ::abort();
}

// ============================================================================
// [mathpresso::mpTraceError]
// ============================================================================

MATHPRESSO_NOAPI Error mpTraceError(Error error) {
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
static void mpDummyFunc(double* result, void*) {
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
struct ContextInternalImpl : public ContextImpl {
  MATHPRESSO_INLINE ContextInternalImpl()
    : _zone(32768 - Zone::kZoneOverhead),
      _heap(&_zone),
      _builder(&_heap),
      _scope(&_builder, static_cast<AstScope*>(NULL), kAstScopeGlobal) {
    mpAtomicSet(&_refCount, 1);
  }
  MATHPRESSO_INLINE ~ContextInternalImpl() {}

  Zone _zone;
  ZoneHeap _heap;
  AstBuilder _builder;
  AstScope _scope;
};

static MATHPRESSO_INLINE ContextImpl* mpContextAddRef(ContextImpl* d) {
  if (d != &mpContextNull)
    mpAtomicInc(&d->_refCount);
  return d;
}

static MATHPRESSO_INLINE void mpContextRelease(ContextImpl* d) {
  if (d != &mpContextNull && !mpAtomicDec(&d->_refCount))
    delete static_cast<ContextInternalImpl*>(d);
}

static ContextImpl* mpContextClone(ContextImpl* otherD_) {
  ContextInternalImpl* d = new(std::nothrow) ContextInternalImpl();
  if (MATHPRESSO_UNLIKELY(d == nullptr)) return nullptr;

  if (otherD_ != &mpContextNull) {
    ContextInternalImpl* otherD = static_cast<ContextInternalImpl*>(otherD_);
    AstSymbolHashIterator it(otherD->_scope._symbols);

    while (it.has()) {
      AstSymbol* sym = it.get();

      StringRef name(sym->_name, sym->_length);
      uint32_t hVal = sym->getHVal();
      uint32_t type = sym->getSymbolType();

      AstSymbol* clonedSym = d->_builder.newSymbol(name, hVal, type, otherD->_scope.getScopeType());
      if (MATHPRESSO_UNLIKELY(clonedSym == nullptr)) {
        delete d;
        return nullptr;
      }

      clonedSym->_symbolFlags = sym->_symbolFlags;
      switch (type) {
        case kAstSymbolVariable:
          clonedSym->setVarSlotId(sym->getVarSlotId());
          clonedSym->setVarOffset(sym->getVarOffset());
          clonedSym->_valueComp = sym->getValueComp();
          break;

        case kAstSymbolIntrinsic:
        case kAstSymbolFunction:
          clonedSym->setOpType(sym->getOpType());
          clonedSym->setFuncArgs(sym->getFuncArgs());
          clonedSym->setFuncPtr(sym->getFuncPtr());
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

static Error mpContextMutable(Context* self, ContextInternalImpl** out) {
  ContextImpl* d = self->_d;

  if (d != &mpContextNull && mpAtomicGet(&d->_refCount) == 1) {
    *out = static_cast<ContextInternalImpl*>(d);
    return kErrorOk;
  }
  else {
    d = mpContextClone(d);
    if (MATHPRESSO_UNLIKELY(d == nullptr))
      return MATHPRESSO_TRACE_ERROR(kErrorNoMemory);

    mpContextRelease(
      mpAtomicSetXchgT<ContextImpl*>(
        &self->_d, d));

    *out = static_cast<ContextInternalImpl*>(d);
    return kErrorOk;
  }
}


// ============================================================================
// [mathpresso::Context - Copy / Reset]
// ============================================================================

Error Context::reset() {
  mpContextRelease(
    mpAtomicSetXchgT<ContextImpl*>(
      &_d, const_cast<ContextImpl*>(&mpContextNull)));

  return kErrorOk;
}

Context& Context::operator=(const Context& other) {
  mpContextRelease(
    mpAtomicSetXchgT<ContextImpl*>(
      &_d, mpContextAddRef(other._d)));
  return *this;
}


// ============================================================================
// [mathpresso::Context - Construction / Destruction]
// ============================================================================

Context::Context()
  : _d(const_cast<ContextImpl*>(&mpContextNull)) {}

Context::Context(const Context& other)
  : _d(mpContextAddRef(other._d)) {}

Context::~Context() {
  mpContextRelease(_d);
}


// ============================================================================
// [mathpresso::Context - Interface]
// ============================================================================

struct GlobalConstant {
  char name[8];
  double value;
};

Error Context::addBuiltIns(void) {
  ContextInternalImpl* d;
  MATHPRESSO_PROPAGATE(mpContextMutable(this, &d));

  uint32_t i;

  for (i = kOpNone + 1; i < kOpCount; i++) 
  {
    const OpInfo& op = OpInfo::get(i);
    MATHPRESSO_ASSERT(op.type == i);

    if (!op.isIntrinsic())
      continue;

    StringRef name(op.name.data());
    uint32_t hVal = HashUtils::hashString(name.getData(), name.getLength());

    AstSymbol* sym = d->_builder.newSymbol(name, hVal, kAstSymbolIntrinsic, kAstScopeGlobal);
    MATHPRESSO_NULLCHECK(sym);

    sym->setDeclared();
    sym->setOpType(op.type);
    sym->setFuncArgs(op.getOpCount());
    sym->setFuncPtr(nullptr);

    d->_scope.putSymbol(sym);
  }

  const GlobalConstant mpGlobalConstants[] = {
    { "NaN", mpGetNan() },
    { "INF", mpGetInf() },
    { "PI" , 3.14159265358979323846 },
    { "E"  , 2.7182818284590452354  }
  };

  for (i = 0; i < MATHPRESSO_ARRAY_SIZE(mpGlobalConstants); i++) {
    const GlobalConstant& c = mpGlobalConstants[i];

    StringRef name(c.name, ::strlen(c.name));
    uint32_t hVal = HashUtils::hashString(name.getData(), name.getLength());

    AstSymbol* sym = d->_builder.newSymbol(name, hVal, kAstSymbolVariable, kAstScopeGlobal);
    MATHPRESSO_NULLCHECK(sym);

    sym->setSymbolFlag(kAstSymbolIsDeclared | kAstSymbolIsAssigned | kAstSymbolIsReadOnly);
    sym->setVarSlotId(kInvalidSlot);
    sym->setVarOffset(0);
    sym->setValue(c.value);

    d->_scope.putSymbol(sym);
  }

  AstSymbol* sym = d->_builder.newSymbol(StringRef("i"), HashUtils::hashString("i",1), kAstSymbolVariable, kAstScopeGlobal);
  MATHPRESSO_NULLCHECK(sym);

  sym->setSymbolFlag(kAstSymbolIsDeclared | kAstSymbolIsAssigned | kAstSymbolIsReadOnly | kAstSymbolIsComplex);
  sym->setVarSlotId(kInvalidSlot);
  sym->setVarOffset(0);
  sym->setValue({0, 1});

  d->_scope.putSymbol(sym);

  return kErrorOk;
}

Error Context::addSymbol(AstSymbol* &sym, const char * name, int type)
{
	ContextInternalImpl* d;
	MATHPRESSO_PROPAGATE(mpContextMutable(this, &d));

	size_t nlen = strlen(name); 
	uint32_t hVal = HashUtils::hashString(name, nlen); 
	sym = d->_scope.getSymbol(StringRef(name, nlen), hVal); 
	if (sym != nullptr) 
		return MATHPRESSO_TRACE_ERROR(kErrorSymbolAlreadyExists); 
    
	sym = d->_builder.newSymbol(StringRef(name, nlen), hVal, type, kAstScopeGlobal);
	if (sym == nullptr) 
		return MATHPRESSO_TRACE_ERROR(kErrorNoMemory); 
	d->_scope.putSymbol(sym); 

	return kErrorOk;
}

Error Context::addConstant(const char* name, double value) 
{
  AstSymbol* sym;
  MATHPRESSO_PROPAGATE(addSymbol(sym, name, kAstSymbolVariable));

  sym->setValue(value);
  sym->setSymbolFlag(kAstSymbolIsDeclared | kAstSymbolIsReadOnly | kAstSymbolIsAssigned);

  return kErrorOk;
}

Error Context::addConstant(const char* name, std::complex<double> value) {
	AstSymbol* sym;
	MATHPRESSO_PROPAGATE(addSymbol(sym, name, kAstSymbolVariable));

	sym->setValue(value);
	sym->setSymbolFlag(kAstSymbolIsDeclared | kAstSymbolIsReadOnly | kAstSymbolIsAssigned | kAstSymbolIsComplex);

	return kErrorOk;
}

Error Context::addVariable(const char* name, int offset, unsigned int flags) 
{
	AstSymbol* sym;
	MATHPRESSO_PROPAGATE(addSymbol(sym, name, kAstSymbolVariable));

	sym->setSymbolFlag(kAstSymbolIsDeclared);
	if (flags & kVariableCplx)
		sym->setSymbolFlag(kAstSymbolIsComplex);

	sym->setVarSlotId(kInvalidSlot);
	sym->setVarOffset(offset);

	if (flags & kVariableRO)
		sym->setSymbolFlag(kAstSymbolIsReadOnly);

	return kErrorOk;
}

Error Context::listSymbols(std::vector<std::string> &syms)
{
	syms.clear();

	ContextInternalImpl* d;
	MATHPRESSO_PROPAGATE(mpContextMutable(this, &d));

	HashIterator<StringRef, AstSymbol> it(d ->_scope.getSymbols());
	do 
	{
		syms.push_back(it.get()->_name);
	} 
	while (it.next());

	return kErrorOk;
}


Error Context::addFunction(const char* name, void* fn, unsigned int flags) 
{
	AstSymbol * sym;
	Error e = addSymbol(sym, name, kAstSymbolFunction);
	if (e != kErrorOk)
	{
		if (e != kErrorSymbolAlreadyExists || sym ->getSymbolType() != kAstSymbolFunction)
			return e;
	}
	else
	{
		sym->setSymbolFlag(kAstSymbolIsDeclared);
	}
  
	// Declaring complex part?
	if (flags & kFunctionTakesComplex)
	{
		if (sym ->getFuncPtr(true))
			return kErrorSymbolAlreadyExists;

		sym ->setFuncPtr(fn, true);
		if (0 == (flags & kFunctionReturnsComplex))
			sym ->setSymbolFlag(kAstSymbolComplexFunctionReturnsReal);
	}
	else
	{
		if (sym ->getFuncPtr())
			return kErrorSymbolAlreadyExists;

		sym ->setFuncPtr(fn);
		if (flags & kFunctionReturnsComplex)
			sym ->setSymbolFlag(kAstSymbolRealFunctionReturnsComplex);
	}

	if (flags & kFunctionHasState)
	  sym->setSymbolFlag(kAstSymbolHasState);

	sym->setFuncArgs(flags & _kFunctionArgMask);
	return kErrorOk;
}

Error Context::delSymbol(const char* name) {
  ContextInternalImpl* d;
  MATHPRESSO_PROPAGATE(mpContextMutable(this, &d));

  size_t nlen = strlen(name);
  uint32_t hVal = HashUtils::hashString(name, nlen);

  AstSymbol* sym = d->_scope.getSymbol(StringRef(name, nlen), hVal);
  if (sym == nullptr)
    return MATHPRESSO_TRACE_ERROR(kErrorSymbolNotFound);

  d->_builder.deleteSymbol(d->_scope.removeSymbol(sym));
  return kErrorOk;
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

Error Expression::compile(const Context& ctx, const char* body, unsigned int options, OutputLog* log) {
	// Init options first.
	options &= _kOptionsMask;

	if (log != nullptr)
		options |= kInternalOptionLog;
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
	{ MATHPRESSO_PROPAGATE(Parser(&ast, &errorReporter, body, len).parseProgram(ast.getProgramNode())); }

	if (options & kOptionDebugAst) {
		ast.dump(sbTmp);
		log->log(OutputLog::kMessageAstInitial, 0, 0, sbTmp.getData(), sbTmp.getLength());
		sbTmp.clear();
	}

	// Perform basic optimizations at AST level.
	{ MATHPRESSO_PROPAGATE(AstOptimizer(&ast, &errorReporter).onProgram(ast.getProgramNode())); }

	if (options & kOptionDebugAst) {
		ast.dump(sbTmp);
		log->log(OutputLog::kMessageAstFinal, 0, 0, sbTmp.getData(), sbTmp.getLength());
		sbTmp.clear();
	}
	
	_isComplex = ast._programNode->returnsComplex();

	// Compile the function to machine code.
	reset();

	CompiledFunc fn = mpCompileFunction(&ast, options, log, _isComplex);

	if (fn == nullptr)
		return MATHPRESSO_TRACE_ERROR(kErrorNoMemory);
	_func = fn;

	return kErrorOk;
}

bool Expression::isCompiled() const 
{
	return _func != mpDummyFunc;
}

void Expression::reset() 
{
	// Allocated by a JIT memory manager, free it.
	if (_func != mpDummyFunc) {
		mpFreeFunction((void*)_func);
		_func = mpDummyFunc;
	}
}

// ============================================================================
// [mathpresso::OutputLog - Construction / Destruction]
// ============================================================================

OutputLog::OutputLog() {}
OutputLog::~OutputLog() {}

// ============================================================================
// [mathpresso::ErrorReporter - Interface]
// ============================================================================

void ErrorReporter::getLineAndColumn(uint32_t position, uint32_t& line, uint32_t& column) {
  // Shouldn't happen, but be defensive.
  if (static_cast<size_t>(position) >= _len) {
    line = 0;
    column = 0;
    return;
  }

  const char* pStart = _body;
  const char* p = pStart + position;

  uint32_t x = 0;
  uint32_t y = 1;

  // Find column.
  while (p[0] != '\n') {
    x++;
    if (p == pStart)
      break;
    p--;
  }

  // Find line.
  while (p != pStart) {
    y += p[0] == '\n';
    p--;
  }

  line = y;
  column = x;
}

void ErrorReporter::onWarning(uint32_t position, const char* fmt, ...) {
  if (reportsWarnings()) {
    StringBuilderTmp<256> sb;

    va_list ap;
    va_start(ap, fmt);

    sb.appendFormatVA(fmt, ap);

    va_end(ap);
    onWarning(position, sb);
  }
}

void ErrorReporter::onWarning(uint32_t position, const StringBuilder& msg) {
  if (reportsWarnings()) {
    uint32_t line, column;
    getLineAndColumn(position, line, column);
    _log->log(OutputLog::kMessageWarning, line, column, msg.getData(), msg.getLength());
  }
}

Error ErrorReporter::onError(Error error, uint32_t position, const char* fmt, ...) {
  if (reportsErrors()) {
    StringBuilderTmp<256> sb;

    va_list ap;
    va_start(ap, fmt);

    sb.appendFormatVA(fmt, ap);

    va_end(ap);
    return onError(error, position, sb);
  }
  else {
    return MATHPRESSO_TRACE_ERROR(error);
  }
}

Error ErrorReporter::onError(Error error, uint32_t position, const StringBuilder& msg) {
  if (reportsErrors()) {
    uint32_t line, column;
    getLineAndColumn(position, line, column);
    _log->log(OutputLog::kMessageError, line, column, msg.getData(), msg.getLength());
  }

  return MATHPRESSO_TRACE_ERROR(error);
}

} // mathpresso namespace
