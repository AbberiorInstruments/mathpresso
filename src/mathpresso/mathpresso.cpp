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
#define ROW(opType, altType, params, precedence, assignment, intrinsic, flags, name) \
  { \
	name, \
    static_cast<uint8_t>(kOp##opType), \
    static_cast<uint8_t>(kOp##altType), \
    static_cast<uint8_t>(precedence), \
    static_cast<uint32_t>( \
      flags | (assignment != 0 ? kOpFlagAssign : 0) \
            | (params     == 1 ? kOpFlagUnary : (params == 2 ? kOpFlagBinary : (params == 3 ? kOpFlagTernary : 0))) \
            | (intrinsic  == 1 ? kOpFlagIntrinsic : 0)), \
  }
#define LTR 0
#define RTL kOpFlagRightToLeft
#define F(flag) kOpFlag##flag
#define CtoC (kOpFlagReturnsComplex | kOpFlagComplex| kOpFlagNoOther)
#define DtoD (kOpFlagNoOther)
#define CtoD (kOpFlagComplex | kOpFlagNoOther)
#define DtoC (kOpFlagReturnsComplex | kOpFlagNoOther)
const OpInfo mpOpInfo[kOpCount] = {
  // +-------------+----------+--+--+--+--+-----+------------------------------------+------------+
  // |Operator     | Alt.     |#N|#P|:=|#I|Assoc| Flags                              | Name       |
  // +-------------+----------+--+--+--+--+-----+------------------------------------+------------+
  OpInfo("<none>", kOpNone, 0, LTR),
  OpInfo("-", kOpNeg, 3, RTL | kOpFlagArithmetic | kOpFlagUnary),
  //ROW(Not          , Not      , 1, 3, 0, 0, RTL | F(Condition)                       , "!"        ),
  OpInfo("!", kOpNot, 3, RTL | kOpFlagCondition | kOpFlagUnary),
  //ROW(IsNan        , IsNan    , 1, 0, 0, 1, LTR | F(Condition)                       , "isnan"    ),
  OpInfo("isnan", kOpIsNan, 0, LTR | kOpFlagCondition | kOpFlagUnary, mpIsNan),
  //ROW(IsInf        , IsInf    , 1, 0, 0, 1, LTR | F(Condition)                       , "isinf"    ),
  OpInfo("isinf", kOpIsInf, 0, LTR | kOpFlagCondition | kOpFlagUnary, mpIsInf),
  //ROW(IsFinite     , IsFinite , 1, 0, 0, 1, LTR | F(Condition)                       , "isfinite" ),
  OpInfo("isfinite", kOpIsFinite, 0, LTR | kOpFlagCondition | kOpFlagUnary, mpIsFinite),
  //ROW(SignBit      , SignBit  , 1, 0, 0, 1, LTR | F(Condition)                       , "signbit"  ),
  OpInfo("signbit", kOpSignBit, 0, LTR | kOpFlagCondition | kOpFlagUnary, mpSignBit),
  //ROW(Round        , Round    , 1, 0, 0, 1, LTR | F(Rounding)                        , "round"    ),
  OpInfo("round", kOpRound, 0, LTR | kOpFlagRounding | kOpFlagUnary, mpRound),
  //ROW(RoundEven    , RoundEven, 1, 0, 0, 1, LTR | F(Rounding)                        , "roundeven"),
  OpInfo("roundeven", kOpRoundEven, 0, LTR | kOpFlagRounding | kOpFlagUnary, mpRoundEven),
  //ROW(Trunc        , Trunc    , 1, 0, 0, 1, LTR | F(Rounding)                        , "trunc"    ),
  OpInfo("trunc", kOpTrunc, 0, LTR | kOpFlagRounding | kOpFlagUnary, mpTrunc),
  //ROW(Floor        , Floor    , 1, 0, 0, 1, LTR | F(Rounding)                        , "floor"    ),
  OpInfo("floor", kOpFloor, 0, LTR | kOpFlagRounding | kOpFlagUnary, mpFloor),
  //ROW(Ceil         , Ceil     , 1, 0, 0, 1, LTR | F(Rounding)                        , "ceil"     ),
  OpInfo("ceil", kOpCeil, 0, LTR | kOpFlagRounding | kOpFlagUnary, mpCeil),
  //ROW(Abs          , Abs      , 1, 0, 0, 1, LTR | 0                                  , "abs"      ),
  OpInfo("abs", kOpAbs, 0, LTR | kOpFlagUnary, mpAbs),
  //ROW(Exp          , Exp      , 1, 0, 0, 1, LTR | 0                                  , "exp"      ),
  OpInfo("exp", kOpExp, 0, LTR | kOpFlagUnary, mpExp, mpExpC),
  //ROW(Log          , Log      , 1, 0, 0, 1, LTR | 0                                  , "log"      ),
  OpInfo("log", kOpLog, 0, LTR | kOpFlagUnary, mpLog, mpLogC),
  //ROW(Log2         , Log2     , 1, 0, 0, 1, LTR | 0                                  , "log2"     ),
  OpInfo("log2", kOpLog2, 0, LTR | kOpFlagUnary, mpLog2, mpLog2C),
  //ROW(Log10        , Log10    , 1, 0, 0, 1, LTR | 0                                  , "log10"    ),
  OpInfo("log10", kOpLog10, 0, LTR | kOpFlagUnary, mpLog10, mpLog10C),
  //ROW(Sqrt         , Sqrt     , 1, 0, 0, 1, LTR | 0                                  , "sqrt"     ),
  OpInfo("sqrt", kOpSqrt, 0, LTR | DtoD | kOpFlagUnary, mpSqrt),
  //ROW(Frac         , Frac     , 1, 0, 0, 1, LTR | 0                                  , "frac"     ),
  OpInfo("frac", kOpFrac, 0, LTR | kOpFlagUnary, mpFrac),
  //ROW(Recip        , Recip    , 1, 0, 0, 1, LTR | 0                                  , "recip"    ),
  OpInfo("recip", kOpRecip, 0, LTR | kOpFlagUnary, mpRecip, mpRecipC),
  //ROW(Sin          , Sin      , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "sin"      ),
  OpInfo("sin", kOpSin, 0, LTR | kOpFlagTrigonometric| kOpFlagUnary, mpSin, mpSinC),
  //ROW(Cos          , Cos      , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "cos"      ),
  OpInfo("cos", kOpCos, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary, mpCos, mpCosC),
  //ROW(Tan          , Tan      , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "tan"      ),
  OpInfo("tan", kOpTan, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary, mpTan, mpTanC),
  //ROW(Sinh         , Sinh     , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "sinh"     ),
  OpInfo("sinh", kOpSinh, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary, mpSinh, mpSinhC),
  //ROW(Cosh         , Cosh     , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "cosh"     ),
  OpInfo("cosh", kOpCosh, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary, mpCosh, mpCoshC),
  //ROW(Tanh         , Tanh     , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "tanh"     ),
  OpInfo("tanh", kOpTanh, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary, mpTanh, mpTanhC),
  //ROW(Asin         , Asin     , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "asin"     ),
  OpInfo("asin", kOpAsin, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary, mpAsin, mpAsinC),
  //ROW(Acos         , Acos     , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "acos"     ),
  OpInfo("acos", kOpAcos, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary, mpAcos, mpAcosC),
  //ROW(Atan         , Atan     , 1, 0, 0, 1, LTR | F(Trigonometric)                   , "atan"     ),
  OpInfo("atan", kOpAtan, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary, mpAtan, mpAtanC),
  //ROW(Assign       , Assign   , 2,15,-1, 0, RTL | 0                                  , "="        ),
  OpInfo("=", kOpAssign, 15, RTL | kOpFlagAssign| kOpFlagBinary),
  //ROW(Eq           , Eq       , 2, 9, 0, 0, LTR | F(Condition)                       , "=="       ),
  OpInfo("==", kOpEq, 9, LTR | kOpFlagCondition | kOpFlagBinary),
  //ROW(Ne           , Ne       , 2, 9, 0, 0, LTR | F(Condition)                       , "!="       ),
  OpInfo("!=", kOpNe, 9, LTR | kOpFlagCondition | kOpFlagBinary),
  //ROW(Lt           , Lt       , 2, 8, 0, 0, LTR | F(Condition)                       , "<"        ),
  OpInfo("<", kOpLt, 8, LTR | kOpFlagCondition | kOpFlagBinary),
  //ROW(Le           , Le       , 2, 8, 0, 0, LTR | F(Condition)                       , "<="       ),
  OpInfo("<=", kOpLe, 8, LTR | kOpFlagCondition | kOpFlagBinary),
  //ROW(Gt           , Gt       , 2, 8, 0, 0, LTR | F(Condition)                       , ">"        ),
  OpInfo(">", kOpGt, 8, LTR | kOpFlagCondition | kOpFlagBinary),
  //ROW(Ge           , Ge       , 2, 8, 0, 0, LTR | F(Condition)                       , ">="       ),
  OpInfo("<=", kOpGe, 8, LTR | kOpFlagCondition | kOpFlagBinary),
  //ROW(Add          , Add      , 2, 6, 0, 0, LTR | F(Arithmetic)    | F(NopIfZero)    , "+"        ),
  OpInfo("+", kOpAdd , 6, LTR | kOpFlagArithmetic | kOpFlagNopIfZero | kOpFlagBinary),
  //ROW(Sub          , Sub      , 2, 6, 0, 0, LTR | F(Arithmetic)    | F(NopIfRZero)   , "-"        ),
  OpInfo("-", kOpSub , 6, LTR | kOpFlagArithmetic | kOpFlagNopIfRZero | kOpFlagBinary),
  //ROW(Mul          , Mul      , 2, 5, 0, 0, LTR | F(Arithmetic)    | F(NopIfOne)     , "*"        ),
  OpInfo("*", kOpMul , 5, LTR | kOpFlagArithmetic | kOpFlagNopIfOne | kOpFlagBinary),
  //ROW(Div          , Div      , 2, 5, 0, 0, LTR | F(Arithmetic)    | F(NopIfROne)    , "/"        ),
  OpInfo("/", kOpDiv , 5, LTR | kOpFlagArithmetic | kOpFlagNopIfROne | kOpFlagBinary),
  //ROW(Mod          , Mod      , 2, 5, 0, 0, LTR | F(Arithmetic)                      , "%"        ),
  OpInfo("%", kOpMod , 5, LTR | DtoD | kOpFlagBinary, mpMod),
  //ROW(Avg          , Avg      , 2, 0, 0, 1, LTR | 0                                  , "avg"      ),
  OpInfo("avg", kOpAvg , 0, LTR | DtoD | kOpFlagBinary, mpAvg),
  //ROW(Min          , Min      , 2, 0, 0, 1, LTR | 0                                  , "min"      ),
  OpInfo("min", kOpMin , 0, LTR | DtoD | kOpFlagBinary, mpMin<double>),
  //ROW(Max          , Max      , 2, 0, 0, 1, LTR | 0                                  , "max"      ),
  OpInfo("max", kOpMax , 0, LTR | DtoD | kOpFlagBinary, mpMax<double>),
  //ROW(Pow          , Pow      , 2, 0, 0, 1, LTR |                    F(NopIfROne)    , "pow"      ),
  OpInfo("pow", kOpPow , 0, LTR | kOpFlagBinary | kOpFlagNopIfROne, mpPow, mpPowC),
  //ROW(Atan2        , Atan2    , 2, 0, 0, 1, LTR | F(Trigonometric)                   , "atan2"    ),
  OpInfo("atan2", kOpAtan2 , 0, LTR | DtoD | kOpFlagTrigonometric | kOpFlagBinary, mpAtan2),
  //ROW(Hypot        , Hypot    , 2, 0, 0, 1, LTR | F(Trigonometric)                   , "hypot"    ),
  OpInfo("hypot", kOpHypot , 0, LTR | DtoD | kOpFlagTrigonometric | kOpFlagBinary, mpHypot),
  //ROW(CopySign     , CopySign , 2, 0, 0, 1, LTR | 0                                  , "copysign" ),
  OpInfo("copysign", kOpCopySign , 0, LTR | DtoD | kOpFlagBinary, mpCopySign),
  //ROW(QMark        , QMark    , 3,15, 0, 0, RTL | 0									 , "?"		  ),
  OpInfo("?", kOpQMark, 15, RTL | kOpFlagTernary),
  //ROW(Colon        , Colon    , 3,15, 0, 0, RTL | 0									 , ":"        ),
  OpInfo(":", kOpColon, 15, RTL | kOpFlagTernary),
  
  //ROW(Real         , Real     , 1, 0, 0, 1, LTR | F(Complex)                         , "real"  ),
  OpInfo("real", kOpReal, 0, LTR | CtoD | kOpFlagUnary , nullptr, mpGetReal),
  //ROW(Imag	       , Imag     , 1, 0, 0, 1, LTR | F(Complex)                         , "imag"  ),
  OpInfo("imag", kOpImag, 0, LTR | CtoD | kOpFlagUnary, nullptr, mpGetImag),
  //ROW(Conjug       , Conjug   , 1, 0, 0, 0, LTR | F(ReturnsComplex) | F(Complex)     , "conjug"   ),
  OpInfo("conjug", kOpConjug, 0, LTR | CtoC | kOpFlagUnary),

  OpInfo("exp_", kOpExpC, 0, LTR | kOpFlagUnary, mpExp, mpExpC),
  //ROW(PowC         , PowC     , 2, 0, 0, 1, LTR | F(ReturnsComplex) | F(Complex)     , "pow_"      ),
  OpInfo("pow_", kOpPowC, 0, LTR | kOpFlagBinary | kOpFlagNopIfROne, mpPow, mpPowC),

  //ROW(LogC         , LogC     , 1, 0, 0, 1, LTR | F(ReturnsComplex) | F(Complex)      , "log_"),
  OpInfo("log_", kOpLogC, 0, LTR | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpLog, mpLogC),
  //ROW(Log2C        , Log2C    , 1, 0, 0, 1, LTR | F(ReturnsComplex) | F(Complex)      , "log2_"),
  OpInfo("log2_", kOpLog2C, 0, LTR | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpLog2, mpLog2C),
  //ROW(Log10C       , Log10C   , 1, 0, 0, 1, LTR | F(ReturnsComplex) | F(Complex)      , "log10_"),
  OpInfo("log10_", kOpLog10C, 0, LTR | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpLog10, mpLog10C),

  //ROW(SqrtC	       , SqrtC    , 1, 0, 0, 1, LTR | F(ReturnsComplex) | F(Complex)     , "sqrtC"    ),
  OpInfo("sqrtC", kOpSqrtC, 0, LTR | CtoC | kOpFlagUnary, nullptr, mpSqrtC),
  //ROW(RecipC       , RecipC   , 1, 0, 0, 1, LTR | F(ReturnsComplex) | F(Complex)     , "recip_"    ),
  OpInfo("recip_", kOpRecipC, 0, LTR | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpRecip, mpRecipC),
  
  //ROW(SinC         , SinC     , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "sin_"),
  OpInfo("sin_", kOpSinC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpSin, mpSinC),
  //ROW(CosC         , CosC     , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "cos_"),
  OpInfo("cos_", kOpCosC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpCos, mpCosC),
  //ROW(TanC         , TanC     , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "tan_"),
  OpInfo("tan_", kOpTanC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpTan, mpTanC),
  
  //ROW(SinhC        , SinhC    , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "sinh_"),
  OpInfo("sinh_", kOpSinhC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpSinh, mpSinhC),
  //ROW(CoshC        , CoshC    , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "cosh_"),
  OpInfo("cosh_", kOpCoshC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpCosh, mpCoshC),
  //ROW(TanhC        , TanhC    , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "tanh_"),
  OpInfo("tanh_", kOpTanhC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpTanh, mpTanhC),

  //ROW(AsinC        , AsinC    , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "asin_"),
  OpInfo("asin_", kOpAsinC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpAsin, mpAsinC),
  //ROW(AcosC        , AcosC    , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "acos_"),
  OpInfo("acos_", kOpAcosC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpAcos, mpAcosC),
  //ROW(AtanC        , AtanC    , 1, 0, 0, 1, LTR | F(Trigonometric) | F(ReturnsComplex) | F(Complex), "atan_")
  OpInfo("atan_", kOpAtanC, 0, LTR | kOpFlagTrigonometric | kOpFlagUnary | kOpFlagReturnsComplex | kOpFlagComplex, mpAtan, mpAtanC)
};
#undef F
#undef RTL
#undef LTR
#undef ROW

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
  if (MATHPRESSO_UNLIKELY(d == NULL)) return NULL;

  if (otherD_ != &mpContextNull) {
    ContextInternalImpl* otherD = static_cast<ContextInternalImpl*>(otherD_);
    AstSymbolHashIterator it(otherD->_scope._symbols);

    while (it.has()) {
      AstSymbol* sym = it.get();

      StringRef name(sym->_name, sym->_length);
      uint32_t hVal = sym->getHVal();
      uint32_t type = sym->getSymbolType();

      AstSymbol* clonedSym = d->_builder.newSymbol(name, hVal, type, otherD->_scope.getScopeType());
      if (MATHPRESSO_UNLIKELY(clonedSym == NULL)) {
        delete d;
        return NULL;
      }

      clonedSym->_symbolFlags = sym->_symbolFlags;
      switch (type) {
        case kAstSymbolVariable:
          clonedSym->setVarSlotId(sym->getVarSlotId());
          clonedSym->setVarOffset(sym->getVarOffset());
          clonedSym->_value = sym->getValue();
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
    if (MATHPRESSO_UNLIKELY(d == NULL))
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
    sym->setFuncPtr(NULL);

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
	if (sym != NULL) 
		return MATHPRESSO_TRACE_ERROR(kErrorSymbolAlreadyExists); 
    
	sym = d->_builder.newSymbol(StringRef(name, nlen), hVal, type, kAstScopeGlobal);
	if (sym == NULL) 
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
  if (sym == NULL)
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

	if (log != NULL)
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
	
	_isComplex = ast._programNode->takesComplex();

	// Compile the function to machine code.
	reset();

	CompiledFunc fn = mpCompileFunction(&ast, options, log, _isComplex);

	if (fn == NULL)
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
