// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include "./mpast_p.h"
#include "./mpcompiler_p.h"
#include "./mpeval_p.h"

namespace mathpresso {
	using namespace asmjit;

	// ============================================================================
	// [mathpresso::JitGlobal]
	// ============================================================================

	struct JitGlobal {
		JitRuntime runtime;
	};
	static JitGlobal jitGlobal;

	// ============================================================================
	// [mathpresso::JitUtils]
	// ============================================================================

	struct JitUtils {
		static void* getFuncByOp(uint32_t op) {
			switch (op) {
			case kOpIsNan: return (void*)(Arg1Func)mpIsNan;
			case kOpIsInf: return (void*)(Arg1Func)mpIsInf;
			case kOpIsFinite: return (void*)(Arg1Func)mpIsFinite;
			case kOpSignBit: return (void*)(Arg1Func)mpSignBit;

			case kOpRound: return (void*)(Arg1Func)mpRound;
			case kOpRoundEven: return (void*)(Arg1Func)rint;
			case kOpTrunc: return (void*)(Arg1Func)trunc;
			case kOpFloor: return (void*)(Arg1Func)floor;
			case kOpCeil: return (void*)(Arg1Func)ceil;

			case kOpAbs: return (void*)(Arg1Func)fabs;
			case kOpExp: return (void*)(Arg1Func)exp;
			case kOpLog: return (void*)(Arg1Func)log;
			case kOpLog2: return (void*)(Arg1Func)log2;
			case kOpLog10: return (void*)(Arg1Func)log10;
			case kOpSqrt: return (void*)(Arg1Func)sqrt;
			case kOpFrac: return (void*)(Arg1Func)mpFrac;
			case kOpRecip: return (void*)(Arg1Func)mpRecip;

			case kOpSin: return (void*)(Arg1Func)sin;
			case kOpCos: return (void*)(Arg1Func)cos;
			case kOpTan: return (void*)(Arg1Func)tan;
			case kOpSinh: return (void*)(Arg1Func)sinh;
			case kOpCosh: return (void*)(Arg1Func)cosh;
			case kOpTanh: return (void*)(Arg1Func)tanh;
			case kOpAsin: return (void*)(Arg1Func)asin;
			case kOpAcos: return (void*)(Arg1Func)acos;
			case kOpAtan: return (void*)(Arg1Func)atan;

			case kOpAvg: return (void*)(Arg2Func)mpAvg;
			case kOpMin: return (void*)(Arg2Func)mpMin<double>;
			case kOpMax: return (void*)(Arg2Func)mpMax<double>;
			case kOpPow: return (void*)(Arg2Func)pow;
			case kOpAtan2: return (void*)(Arg2Func)atan2;
			case kOpHypot: return (void*)(Arg2Func)hypot;
			case kOpCopySign: return (void*)(Arg2Func)mpCopySign;
			case kOpReal: return (void*)(ArgFuncCtoD)mpGetReal;
			case kOpImag: return (void*)(ArgFuncCtoD)mpGetImag;

			case kOpExpC: return (void*)(ArgFuncC)mpExpC;
			case kOpPowC: return (void*)(ArgFuncC)mpPowC;

			case kOpLogC: return (void*)(ArgFuncC)mpLogC;
			case kOpLog2C: return (void*)(ArgFuncC)mpLog2C;
			case kOpLog10C: return (void*)(ArgFuncC)mpLog10C;

			case kOpSqrtC: return (void*)(ArgFuncC)mpSqrtC;
			case kOpRecipC: return (void*)(ArgFuncC)mpRecipC;

			case kOpSinC: return (void*)(ArgFuncC)mpSinC;
			case kOpCosC: return (void*)(ArgFuncC)mpCosC;
			case kOpTanC: return (void*)(ArgFuncC)mpTanC;

			case kOpSinhC: return (void*)(ArgFuncC)mpSinhC;
			case kOpCoshC: return (void*)(ArgFuncC)mpCoshC;
			case kOpTanhC: return (void*)(ArgFuncC)mpTanhC;

			case kOpAsinC: return (void*)(ArgFuncC)mpAsinC;
			case kOpAcosC: return (void*)(ArgFuncC)mpAcosC;
			case kOpAtanC: return (void*)(ArgFuncC)mpAtanC;

			default:
				MATHPRESSO_ASSERT_NOT_REACHED();
				return NULL;
			}
		}
	};

	// ============================================================================
	// [mathpresso::JitVar]
	// ============================================================================

	struct MATHPRESSO_NOAPI JitVar {
		enum FLAGS {
			FLAG_NONE = 0,
			FLAG_RO = 1
		};

		MATHPRESSO_INLINE JitVar() : op(), flags(FLAG_NONE) {}
		MATHPRESSO_INLINE JitVar(Operand op, uint32_t flags) : op(op), flags(flags) {}
		MATHPRESSO_INLINE JitVar(const JitVar& other) : op(other.op), flags(other.flags) {}
		MATHPRESSO_INLINE ~JitVar() {}

		// Reset
		MATHPRESSO_INLINE void reset() {
			op.reset();
			flags = FLAG_NONE;
		}

		// Operator Overload.
		MATHPRESSO_INLINE const JitVar& operator=(const JitVar& other) {
			op = other.op;
			flags = other.flags;
			return *this;
		}

		// Swap.
		MATHPRESSO_INLINE void swapWith(JitVar& other) {
			JitVar t(*this);
			*this = other;
			other = t;
		}

		// Operand management.
		MATHPRESSO_INLINE const Operand& getOperand() const { return op; }
		MATHPRESSO_INLINE const X86Mem& getMem() const { return *static_cast<const X86Mem*>(&op); }
		MATHPRESSO_INLINE const X86Xmm& getXmm() const { return *static_cast<const X86Xmm*>(&op); }

		MATHPRESSO_INLINE bool isNone() const { return op.isNone(); }
		MATHPRESSO_INLINE bool isMem() const { return op.isMem(); }
		MATHPRESSO_INLINE bool isXmm() const { return op.isReg(X86Reg::kRegXmm); }

		// Flags.
		MATHPRESSO_INLINE bool isRO() const { return (flags & FLAG_RO) != 0; }
		MATHPRESSO_INLINE void setRO() { flags |= FLAG_RO; }
		MATHPRESSO_INLINE void clearRO() { flags &= ~FLAG_RO; }

		// Members.
		Operand op;
		uint32_t flags;
	};

	// ============================================================================
	// [mathpresso::JitCompiler]
	// ============================================================================

	struct MATHPRESSO_NOAPI JitCompiler {
		JitCompiler(ZoneHeap* heap, X86Compiler* c);
		~JitCompiler();

		// Function Generator.
		void beginFunction();
		void endFunction();

		// Function generator - complex.
		void beginFunctionComp();

		// Variable Management.
		JitVar copyVar(const JitVar& other, uint32_t flags);
		JitVar writableVar(const JitVar& other);
		JitVar registerVar(const JitVar& other);

		JitVar copyVarComplex(const JitVar & other, uint32_t flags);
		JitVar writableVarComplex(const JitVar & other);
		JitVar registerVarComplex(const JitVar & other, bool otherIsComplex = false);
		JitVar registerVarAsComplex(const JitVar & other);

		// Compiler.
		void compile(AstBlock* node, AstScope* rootScope, uint32_t numSlots, bool b_complex);
		
		JitVar onNode(AstNode* node);
		JitVar onBlock(AstBlock* node);
		JitVar onVarDecl(AstVarDecl* node);
		JitVar onVar(AstVar* node);
		JitVar onImm(AstImm* node);
		JitVar onUnaryOp(AstUnaryOp* node);
		JitVar onBinaryOp(AstBinaryOp* node);
		JitVar onTernaryOp(AstTernaryOp * node);
		JitVar onCall(AstCall* node);
		JitVar onCallComp(AstCall * node);


		// Helpers.
		void inlineRound(const X86Xmm& dst, const X86Xmm& src, uint32_t op, bool isComplex, bool returnsComplex);
		void inlineCallAbstract(const X86Xmm& dst, const X86Xmm* args, uint32_t count, uint32_t op, bool paramsAreComplex, bool returnsComplex);
		void inlineCall(const X86Xmm& dst, const X86Xmm* args, uint32_t count, void* fn);
		void inlineCallDRetC(const X86Xmm & dst, const X86Xmm * args, uint32_t count, void * fn);
		void inlineCallCRetD(const X86Xmm & dst, const X86Xmm * args, const uint32_t count, void * fn);
		void inlineCallComplex(const X86Xmm & dst, const X86Xmm * args, uint32_t count, void * fn);

		// Constants.
		void prepareConstPool();
		void prepareConstPoolComp();
		JitVar getConstantU64(uint64_t value);
		JitVar getConstantU64Compl(uint64_t real, uint64_t imag);
		JitVar getConstantU64AsPD(uint64_t value);
		JitVar getConstantD64(double value);
		JitVar getConstantD64Compl(std::complex<double> value);
		JitVar getConstantD64AsPD(double value);

		// Members.
		ZoneHeap* heap;
		X86Compiler* cc;

		X86Gp resultAddress;
		X86Gp variablesAddress;

		JitVar* varSlots;
		CBNode* functionBody;

		Label constLabel;
		X86Gp constPtr;
		ConstPool constPool;

		bool enableSSE4_1;
		bool complex = false;
	};

	JitCompiler::JitCompiler(ZoneHeap* heap, X86Compiler* cc)
		: heap(heap),
		cc(cc),
		varSlots(NULL),
		functionBody(NULL),
		constPool(&cc->_cbDataZone) {

		enableSSE4_1 = CpuInfo::getHost().hasFeature(CpuInfo::kX86FeatureSSE4_1);
	}
	JitCompiler::~JitCompiler() {}

	void JitCompiler::beginFunction() {
		cc->addFunc(FuncSignature2<void, double*, double*>(CallConv::kIdHostCDecl));

		resultAddress = cc->newIntPtr("pResult");
		variablesAddress = cc->newIntPtr("pVariables");
		constPtr = cc->newIntPtr("pConst");

		cc->setArg(0, resultAddress);
		cc->setArg(1, variablesAddress);

		functionBody = cc->getCursor();
	}


	//! not necessary, but for completeness.
	void JitCompiler::beginFunctionComp() {
		cc->addFunc(FuncSignature2<void, double*, double*>(CallConv::kIdHostCDecl));

		resultAddress = cc->newIntPtr("pResult");
		variablesAddress = cc->newIntPtr("pVariables");
		constPtr = cc->newIntPtr("pConst");

		cc->setArg(0, resultAddress);
		cc->setArg(1, variablesAddress);

		functionBody = cc->getCursor();
	}

	void JitCompiler::endFunction() {
		cc->endFunc();

		if (constLabel.isValid())
			cc->embedConstPool(constLabel, constPool);
	}

	JitVar JitCompiler::copyVar(const JitVar& other, uint32_t flags) {
		JitVar v(cc->newXmmSd(), flags);
		cc->emit(other.isXmm() ? X86Inst::kIdMovapd : X86Inst::kIdMovsd,
			v.getXmm(), other.getOperand());
		return v;
	}

	JitVar JitCompiler::writableVar(const JitVar& other) {
		if (other.isMem() || other.isRO())
			return copyVar(other, other.flags & ~JitVar::FLAG_RO);
		else
			return other;
	}

	JitVar JitCompiler::registerVar(const JitVar& other) {
		if (other.isMem())
			return copyVar(other, other.flags);
		else
			return other;
	}

	JitVar JitCompiler::copyVarComplex(const JitVar& other, uint32_t flags) {
		JitVar v(cc->newXmmPd(), flags);

		cc->emit(X86Inst::kIdMovupd, v.getXmm(), other.getOperand());
		return v;
	}

	JitVar JitCompiler::writableVarComplex(const JitVar& other) {
		if (other.isMem() || other.isRO())
			return copyVarComplex(other, other.flags & ~JitVar::FLAG_RO);
		else
			return other;
	}

	JitVar JitCompiler::registerVarComplex(const JitVar& other, bool otherIsNonComplex) {
		if (otherIsNonComplex) 
			return registerVarAsComplex(other);
		else if (other.isMem())
			return copyVarComplex(other, other.flags);
		else
			return other;
	}

	JitVar JitCompiler::registerVarAsComplex(const JitVar& other) {
		JitVar v(cc->newXmmPd(), other.flags);
		cc->pxor(v.getXmm(), v.getXmm());
		if (other.isMem()) {
			cc->movlpd(v.getXmm(), other.getMem());
		}
		else {
			cc->movsd(v.getXmm(), other.getXmm());
		}
		return v;
	}

	//! Compiles an AstBlock into assembler.
	//! NOTE: use beginFunction() before and endFunction() after calling this.
	void JitCompiler::compile(AstBlock* node, AstScope* rootScope, uint32_t numSlots, bool b_complex) {
		// Create Definitions for the Variables and add them as JitVar
		if (numSlots != 0) {
			varSlots = static_cast<JitVar*>(heap->alloc(sizeof(JitVar) * numSlots));
			if (varSlots == NULL) return;

			for (uint32_t i = 0; i < numSlots; i++)
				varSlots[i] = JitVar();
		}

		// Result of the function or NaN. Here the AST is compiled.
		JitVar result = onBlock(node);

		// Write altered global variables.
		{
			AstSymbolHashIterator it(rootScope->getSymbols());
			while (it.has()) {
				AstSymbol* sym = it.get();
				if (sym->isGlobal() && sym->isAltered()) {
					JitVar v = varSlots[sym->getVarSlotId()];
					if (!b_complex)
					{
						cc->emit(X86Inst::kIdMovsd,
							x86::ptr(variablesAddress, sym->getVarOffset()), registerVar(v).getXmm());
					}
					else
					{
						cc->emit(X86Inst::kIdMovapd,
							x86::ptr(variablesAddress, sym->getVarOffset()), registerVarComplex(v).getXmm());
					}
				}

				it.next();
			}
		}

		// Return NaN if no result is given.
		X86Xmm var;
		if (!b_complex)
		{
			if (result.isNone())
				var = registerVar(getConstantD64(mpGetNan())).getXmm();
			else
				var = registerVar(result).getXmm();
			cc->movsd(x86::ptr(resultAddress), var);
		}
		else
		{
			if (result.isNone())
				var = registerVarComplex(getConstantD64Compl(mpGetNan())).getXmm();
			else
				var = registerVarComplex(result, !node->hasNodeFlag(kAstComplex)).getXmm();
			cc->movupd(x86::ptr(resultAddress), var);
		}

		// Release the Space allocated for the variables
		if (numSlots != 0)
			heap->release(varSlots, sizeof(JitVar) * numSlots);
	}

	JitVar JitCompiler::onNode(AstNode* node) {
		switch (node->getNodeType()) {
		case kAstNodeBlock: return onBlock(static_cast<AstBlock*>(node));
		case kAstNodeVarDecl: return onVarDecl(static_cast<AstVarDecl*>(node));
		case kAstNodeVar: return onVar(static_cast<AstVar*>(node));
		case kAstNodeImm: return onImm(static_cast<AstImm*>(node));
		case kAstNodeUnaryOp: return onUnaryOp(static_cast<AstUnaryOp*>(node));
		case kAstNodeBinaryOp: return onBinaryOp(static_cast<AstBinaryOp*>(node));
		case kAstNodeTernaryOp: return onTernaryOp(static_cast<AstTernaryOp*>(node));
		case kAstNodeCall: return onCall(static_cast<AstCall*>(node));

		default:
			MATHPRESSO_ASSERT_NOT_REACHED();
			return JitVar();
		}
	}

	JitVar JitCompiler::onBlock(AstBlock* node) {
		JitVar result;
		uint32_t i, len = node->getLength();

		for (i = 0; i < len; i++)
			result = onNode(node->getAt(i));

		// Return the last result (or no result if the block is empty).
		return result;
	}

	JitVar JitCompiler::onVarDecl(AstVarDecl* node) {
		JitVar result;

		if (node->hasChild())
			result = onNode(node->getChild());

		AstSymbol* sym = node->getSymbol();
		uint32_t slotId = sym->getVarSlotId();

		result.setRO();
		varSlots[slotId] = result;

		return result;
	}

	JitVar JitCompiler::onVar(AstVar* node) {
		AstSymbol* sym = node->getSymbol();
		uint32_t slotId = sym->getVarSlotId();
		bool b_complex = node ->hasNodeFlag(kAstReturnsComplex);

		JitVar result = varSlots[slotId];
		if (result.isNone()) {
			if (sym->isGlobal()) {
				result = JitVar(x86::ptr(variablesAddress, sym->getVarOffset()), JitVar::FLAG_RO);
				varSlots[slotId] = result;
				if (sym->getWriteCount() > 0)
				{
					if (!b_complex)
						result = copyVar(result, JitVar::FLAG_NONE);
					else
						result = copyVarComplex(result, JitVar::FLAG_NONE);
				}

			}
			else 
			{
				if (!b_complex)
				{
					result = getConstantD64(mpGetNan());
				}
				else
				{
					result = getConstantD64Compl(mpGetNan());
				}
				varSlots[slotId] = result;
			}
		}

		return result;
	}
	
	JitVar JitCompiler::onImm(AstImm* node) {
		if (node ->isComplex())
			return getConstantD64Compl(node->getValue());
		else
			return getConstantD64(node->getValue());
	}

	JitVar JitCompiler::onUnaryOp(AstUnaryOp* node) {
		uint32_t op = node->getOp();
		JitVar var = onNode(node->getChild());

		if (!node->hasNodeFlag(kAstComplex))
		{
			switch (op)
			{
			case kOpNone:
				return var;

			case kOpNeg:
			{
				X86Xmm result = cc->newXmmSd();
				cc->emit(X86Inst::kIdXorpd, result, result);
				cc->emit(X86Inst::kIdSubsd, result, var.getOperand());
				return JitVar(result, JitVar::FLAG_NONE);
			}

			case kOpNot:
			{
				var = writableVar(var);
				cc->cmpsd(var.getXmm(), getConstantD64AsPD(0.0).getMem(), int(x86::kCmpEQ));
				cc->andpd(var.getXmm(), getConstantD64AsPD(1.0).getMem());
				return var;
			}

			case kOpIsNan:
			{
				var = writableVar(var);
				cc->cmpsd(var.getXmm(), var.getXmm(), int(x86::kCmpEQ));
				cc->andnpd(var.getXmm(), getConstantD64AsPD(1.0).getMem());
				return var;
			}

			case kOpIsInf:
			{
				var = writableVar(var);
				cc->orpd(var.getXmm(), getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
				cc->cmpsd(var.getXmm(), getConstantU64(MATHPRESSO_UINT64_C(0xFF80000000000000)).getMem(), int(x86::kCmpEQ));
				cc->andpd(var.getXmm(), getConstantD64AsPD(1.0).getMem());
				return var;
			}

			case kOpIsFinite:
			{
				var = writableVar(var);
				cc->orpd(var.getXmm(), getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
				cc->cmpsd(var.getXmm(), getConstantD64(0.0).getMem(), int(x86::kCmpLE));
				cc->andpd(var.getXmm(), getConstantD64AsPD(1.0).getMem());
				return var;
			}

			case kOpSignBit:
			{
				X86Xmm result = cc->newXmmSd();
				cc->pshufd(result, registerVar(var).getXmm(), x86::shufImm(3, 2, 1, 1));
				cc->psrad(result, 31);
				cc->andpd(result, getConstantD64AsPD(1.0).getMem());
				return JitVar(result, JitVar::FLAG_NONE);
			}

			case kOpRound:
			case kOpRoundEven:
			case kOpTrunc:
			case kOpFloor:
			case kOpCeil:
			{
				var = writableVar(var);
				inlineRound(var.getXmm(), var.getXmm(), op, node->hasNodeFlag(kAstComplex), node->getChild()->hasNodeFlag(kAstReturnsComplex));
				return var;
			}

			case kOpAbs:
			{
				X86Xmm result = cc->newXmmSd();
				cc->emit(X86Inst::kIdXorpd, result, result);
				cc->emit(X86Inst::kIdSubsd, result, var.getOperand());
				cc->emit(X86Inst::kIdMaxsd, result, var.getOperand());
				return JitVar(result, JitVar::FLAG_NONE);
			}

			case kOpExp:
			case kOpLog:
			case kOpLog2:
			case kOpLog10:
				break;

			case kOpSqrt:
			{
				X86Xmm result = cc->newXmmSd();
				cc->emit(X86Inst::kIdSqrtsd, result, var.getOperand());
				return JitVar(result, JitVar::FLAG_NONE);
			}

			case kOpFrac:
			{
				var = writableVar(var);
				X86Xmm tmp = cc->newXmmSd();

				if (enableSSE4_1)
				{
					cc->emit(X86Inst::kIdRoundsd, tmp, var.getOperand(), int(x86::kRoundDown | x86::kRoundInexact));
					cc->emit(X86Inst::kIdSubsd, var.getOperand(), tmp);
					return var;
				}
				else
				{
					// Pure SSE2 `frac()`, uses the same rounding trick as `floor()`.
					inlineRound(tmp, var.getXmm(), kOpFloor, node->hasNodeFlag(kAstComplex), node->hasNodeFlag(kAstReturnsComplex));
					cc->subsd(var.getXmm(), tmp);
					return var;
				}
			}

			case kOpRecip:
			{
				X86Xmm result = cc->newXmmSd();
				cc->emit(X86Inst::kIdMovsd, result, getConstantD64(1.0).getOperand());
				cc->emit(X86Inst::kIdDivsd, result, var.getOperand());
				return JitVar(result, JitVar::FLAG_NONE);
			}

			case kOpSin:
			case kOpCos:
			case kOpTan:
			case kOpSinh:
			case kOpCosh:
			case kOpTanh:
			case kOpAsin:
			case kOpAcos:
			case kOpAtan:
				break;
			}


			// No inline implementation -> function call.
			X86Xmm result;
			X86Xmm args[1];
			if (node->hasNodeFlag(kAstComplex)) {
				result = cc->newXmmPd();
				args[0] = registerVarComplex(var).getXmm();
			}
			else {
				result = cc->newXmmSd();
				args[0] = registerVar(var).getXmm();
			}


			inlineCallAbstract(result, args, 1, op, node->hasNodeFlag(kAstComplex), node->hasNodeFlag(kAstReturnsComplex));
			return JitVar(result, JitVar::FLAG_NONE);
		
		}
		else {
			switch (op)
			{
			case kOpConjug:
				JitVar result = registerVarComplex(var, node->getChild()->hasNodeFlag(kAstReturnsComplex));
				cc->pxor(result.getXmm(), getConstantU64Compl(uint64_t(0x0000000000000000), uint64_t(0x8000000000000000)).getMem());
				return result;
			}

			X86Xmm result;
			X86Xmm args[1];

			if (node->hasNodeFlag(kAstComplex)) {
				result = cc->newXmmPd();
				args[0] = registerVarComplex(var, !node->getChild()->hasNodeFlag(kAstReturnsComplex)).getXmm();
			}
			else {
				result = cc->newXmmSd();
				args[0] = registerVar(var).getXmm();
			}


			inlineCallAbstract(result, args, 1, op, node->hasNodeFlag(kAstComplex), node->hasNodeFlag(kAstReturnsComplex));
			return JitVar(result, JitVar::FLAG_NONE);
		}
	}

	JitVar JitCompiler::onBinaryOp(AstBinaryOp* node) {
		uint32_t op = node->getOp();

		AstNode* left = node->getLeft();
		AstNode* right = node->getRight();



		// Compile assignment. Should never be reached.
		if (op == kOpAssign) {
			AstVar* varNode = reinterpret_cast<AstVar*>(left);
			MATHPRESSO_ASSERT(varNode->getNodeType() == kAstNodeVar);

			if (varNode->isComplex()) 
			{
				;
			}
			else 
			{
				AstSymbol* sym = varNode->getSymbol();
				uint32_t slotId = sym->getVarSlotId();

				JitVar result = onNode(right);
				result.setRO();

				sym->setAltered();
				varSlots[slotId] = result;
				return result;
			}
		}


		JitVar vl, vr;

		if (left->getNodeType() == kAstNodeVar &&
			right->getNodeType() == kAstNodeVar &&
			static_cast<AstVar*>(left)->getSymbol() == static_cast<AstVar*>(right)->getSymbol())
		{
			vl = vr = writableVar(onNode(node->getLeft()));
		}
		else
		{
			vl = onNode(node->getLeft());
			vr = onNode(node->getRight());

			// Commutativity.
			if (op == kOpAdd || op == kOpMul || op == kOpAvg || op == kOpMin || op == kOpMax)
			{
				if (vl.isRO() && !vr.isRO())
					vl.swapWith(vr);
			}
		
		}

		uint32_t inst = 0;
		int predicate = 0;

		if (!node->hasNodeFlag(kAstComplex))
		{
			vl = writableVar(vl);
			switch (op)
			{
			case kOpEq: predicate = x86::kCmpEQ; goto emitCompare;
			case kOpNe: predicate = x86::kCmpNEQ; goto emitCompare;
			case kOpGt: predicate = x86::kCmpNLE; goto emitCompare;
			case kOpGe: predicate = x86::kCmpNLT; goto emitCompare;
			case kOpLt: predicate = x86::kCmpLT; goto emitCompare;
			case kOpLe: predicate = x86::kCmpLE; goto emitCompare;
			emitCompare:
			{
				vl = writableVar(vl);
				cc->emit(X86Inst::kIdCmpsd, vl.getXmm(), vr.getOperand(), predicate);
				cc->emit(X86Inst::kIdAndpd, vl.getXmm(), getConstantD64AsPD(1.0).getOperand());
				return vl;
			}

			case kOpAdd: inst = X86Inst::kIdAddsd; goto emitInst;
			case kOpSub: inst = X86Inst::kIdSubsd; goto emitInst;
			case kOpMul: inst = X86Inst::kIdMulsd; goto emitInst;
			case kOpDiv: inst = X86Inst::kIdDivsd; goto emitInst;
			case kOpMin: inst = X86Inst::kIdMinsd; goto emitInst;
			case kOpMax: inst = X86Inst::kIdMaxsd; goto emitInst;
			emitInst:
			{
				cc->emit(inst, vl.getOperand(), vr.getOperand()); return vl;
			}

			case kOpAvg:
			{
				vl = writableVar(vl);
				cc->emit(X86Inst::kIdAddsd, vl.getXmm(), vr.getOperand());
				cc->emit(X86Inst::kIdMulsd, vl.getXmm(), getConstantD64(0.5).getOperand());
				return vl;
			}

			case kOpMod:
			{
				X86Xmm result = cc->newXmmSd();
				X86Xmm tmp = cc->newXmmSd();

				vl = writableVar(vl);
				vr = registerVar(vr);

				cc->movsd(result, vl.getXmm());
				cc->divsd(vl.getXmm(), vr.getXmm());
				inlineRound(vl.getXmm(), vl.getXmm(), kOpTrunc, node->hasNodeFlag(kAstComplex), node->hasNodeFlag(kAstReturnsComplex));
				cc->mulsd(vl.getXmm(), vr.getXmm());
				cc->subsd(result, vl.getXmm());

				return JitVar(result, JitVar::FLAG_NONE);
			}

			case kOpPow:
			case kOpAtan2:
			case kOpHypot:
				break;

			case kOpCopySign:
			{
				vl = writableVar(vl);
				vr = writableVar(vr);

				cc->andpd(vl.getXmm(), getConstantU64AsPD(MATHPRESSO_UINT64_C(0x7FFFFFFFFFFFFFFF)).getMem());
				cc->andpd(vr.getXmm(), getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
				cc->orpd(vl.getXmm(), vr.getXmm());

				return vl;
			}

			default:
				MATHPRESSO_ASSERT_NOT_REACHED();
				return vl;
			}
		}
		else
		{

			// Handle the case that the operands are the same variable.
			if (left->isComplex() && right->isComplex() &&
				static_cast<AstVar*>(left)->getSymbol() == static_cast<AstVar*>(right)->getSymbol())
			{
				vl = vr = writableVarComplex(onNode(left));
			}
			else
			{

				if (!left->hasNodeFlag(kAstReturnsComplex)&& vl.isMem())
					vl = registerVarAsComplex(vl);
				else 
					vl = writableVarComplex(vl);

				if (!right->hasNodeFlag(kAstReturnsComplex) && vr.isMem())
					vr = registerVarAsComplex(vr);

				if (vl.isRO() && !vr.isRO() && op == kOpAdd && op == kOpMul)
					vl.swapWith(vr);
				

			}

			if (vr.getOperand() == vl.getOperand())
			{
				vr = copyVarComplex(vl, JitVar::FLAG_NONE);
			}
			else if (vr.getOperand().isMem())
			{
				vr = copyVarComplex(vr, JitVar::FLAG_NONE);
			}


			JitVar ret(cc->newXmmPd(), JitVar::FLAG_NONE);
			JitVar negateImag = getConstantU64Compl(uint64_t(0x0000000000000000), uint64_t(0x8000000000000000));
			
			switch (op) 
			{
			case kOpAdd:
				cc->addpd(vl.getXmm(), vr.getXmm());
				return vl;
			case kOpSub:
				cc->subpd(vl.getXmm(), vr.getXmm());
				return vl;
			case kOpMul:
				cc->movapd(ret.getXmm(), vl.getXmm());
				cc->mulpd(ret.getXmm(), vr.getXmm());
				cc->shufpd(vr.getXmm(), vr.getXmm(), 1);
				cc->pxor(vr.getXmm(), negateImag.getMem());
				cc->mulpd(vl.getXmm(), vr.getXmm());
				cc->hsubpd(ret.getXmm(), vl.getXmm());
				return ret;
			case kOpDiv:
				cc->pxor(vr.getXmm(), negateImag.getMem());

				cc->movapd(ret.getXmm(), vl.getXmm());
				cc->mulpd(ret.getXmm(), vr.getXmm());
				cc->shufpd(vr.getXmm(), vr.getXmm(), 1);
				cc->pxor(vr.getXmm(), negateImag.getMem());
				cc->mulpd(vl.getXmm(), vr.getXmm());
				cc->hsubpd(ret.getXmm(), vl.getXmm());

				cc->mulpd(vr.getXmm(), vr.getXmm());
				cc->haddpd(vr.getXmm(), vr.getXmm());
				cc->divpd(ret.getXmm(), vr.getXmm());
				return ret;
			case kOpPowC:
				break;

			default:
				MATHPRESSO_ASSERT_NOT_REACHED();
				return vl;
			}
		}

		// No inline implementation -> function call.
		X86Xmm result;
		X86Xmm args[2];
		if (node->hasNodeFlag(kAstComplex))
		{
			result = cc->newXmmPd();
			args[0] = registerVarComplex(vl).getXmm();
			args[1] = registerVarComplex(vr).getXmm();
		}
		else
		{
			result = cc->newXmmSd();
			args[0] = registerVar(vl).getXmm();
			args[1] = registerVar(vr).getXmm();
		}

		inlineCallAbstract(result, args, 2, op, node->hasNodeFlag(kAstComplex), node->hasNodeFlag(kAstReturnsComplex));

		return JitVar(result, JitVar::FLAG_NONE);
	}



	JitVar JitCompiler::onTernaryOp(AstTernaryOp* node)
	{
		Label lblElse = cc->newLabel();
		Label lblEnd = cc->newLabel();
		JitVar erg;
		AstNode* left = node->getLeft();
		AstNode* right = node->getRight();
		AstNode* condition = node->getCondition();

		JitVar ret = onNode(condition);
					
		if (condition->hasNodeFlag(kAstReturnsComplex))
			cc->haddpd(ret.getXmm(), ret.getXmm());
			
		cc->ucomisd(ret.getXmm(), getConstantD64(0).getMem());
		cc->je(lblElse);

		X86Xmm regErg = cc->newXmmPd();
		JitVar ergLeft = onNode(left);

		bool lIsVarOrImm = left  ->getNodeType() == kAstNodeVar || left  ->getNodeType() == kAstNodeImm;
		bool rIsVarOrImm = right ->getNodeType() == kAstNodeVar || right ->getNodeType() == kAstNodeImm;

		if (lIsVarOrImm)
		{
			if (left ->isComplex())
			{
				cc->movupd(regErg, ergLeft.getXmm());
			}
			else
			{
				cc->xorpd(regErg, regErg);
				cc->movsd(regErg, ergLeft.getXmm());
			}
		}

		cc->jmp(lblEnd);
		cc->bind(lblElse);

		JitVar ergRight = onNode(right);

		if (rIsVarOrImm)
		{
			if (right ->isComplex())
			{
				cc->movupd(regErg, ergRight.getXmm());
			}
			else
			{
				cc->xorpd(regErg, regErg);
				cc->movsd(regErg, ergRight.getXmm());
			}
		}

		cc->bind(lblEnd);

		return copyVarComplex(JitVar(regErg, JitVar::FLAG_NONE), JitVar::FLAG_NONE);
	}

	JitVar JitCompiler::onCall(AstCall* node)
	{
		if (node->hasNodeFlag(kAstComplex))
		{
			return onCallComp(node);
		}
		uint32_t i, count = node->getLength();
		AstSymbol* sym = node->getSymbol();

		X86Xmm result = node->hasNodeFlag(kAstReturnsComplex) ? cc->newXmmPd() : cc->newXmmSd();
		X86Xmm args[8];

		for (i = 0; i < count; i++)
		{
			args[i] = registerVar(onNode(node->getAt(i))).getXmm();
		}

		if (node->hasNodeFlag(kAstReturnsComplex))
			inlineCallDRetC(result, args, count, sym->getFuncPtr());
		else
			inlineCall(result, args, count, sym->getFuncPtr());
		return JitVar(result, JitVar::FLAG_NONE);
	}

	JitVar JitCompiler::onCallComp(AstCall* node)
	{
		uint32_t i, count = node->getLength();
		AstSymbol* sym = node->getSymbol();

		X86Xmm result = node->hasNodeFlag(kAstReturnsComplex) ? cc->newXmmPd() : cc->newXmmSd();
		X86Xmm args[8];

		for (i = 0; i < count; i++) {
			args[i] = registerVarComplex(onNode(node->getAt(i)), !node->getAt(i)->hasNodeFlag(kAstReturnsComplex)).getXmm();
		}
		if (node->hasNodeFlag(kAstReturnsComplex))
			inlineCallComplex(result, args, count, sym->getFuncPtr(true));
		else
			inlineCallCRetD(result, args, count, sym->getFuncPtr(true));

		return JitVar(result, JitVar::FLAG_NONE);
	}

	void JitCompiler::inlineRound(const X86Xmm& dst, const X86Xmm& src, uint32_t op, bool isComplex, bool returnsComplex) {
		// SSE4.1 implementation is easy except `round()`, which is not `roundeven()`.
		if (enableSSE4_1) {
			if (op == kOpRound) {
				X86Xmm tmp = cc->newXmmSd();
				cc->roundsd(tmp, src, x86::kRoundDown | x86::kRoundInexact);

				if (dst.getId() != src.getId())
					cc->movsd(dst, src);

				cc->subsd(dst, tmp);
				cc->cmpsd(dst, getConstantD64(0.5).getMem(), x86::kCmpNLT);
				cc->andpd(dst, getConstantD64AsPD(1.0).getMem());
				cc->addpd(dst, tmp);
			}
			else {
				int predicate = 0;

				switch (op) {
				case kOpRoundEven: predicate = x86::kRoundNearest; break;
				case kOpTrunc: predicate = x86::kRoundTrunc; break;
				case kOpFloor: predicate = x86::kRoundDown; break;
				case kOpCeil: predicate = x86::kRoundUp; break;
				default:
					MATHPRESSO_ASSERT_NOT_REACHED();
				}
				cc->roundsd(dst, src, predicate | x86::kRoundInexact);
			}
			return;
		}

		// Pure SSE2 requires the following rounding trick:
		//   double roundeven(double x) {
		//     double maxn   = pow(2, 52);
		//     double magic0 = pow(2, 52) + pow(2, 51);
		//     return x >= maxn ? x : x + magic0 - magic0;
		//   }
		const double maxn = 4503599627370496.0;
		const double magic0 = 6755399441055744.0;
		const double magic1 = 6755399441055745.0;

		if (op == kOpRoundEven) {
			X86Xmm t1 = cc->newXmmSd();
			X86Xmm t2 = cc->newXmmSd();

			cc->movsd(t1, src);
			cc->movsd(t2, src);

			cc->addsd(t1, getConstantD64(magic0).getMem());
			cc->cmpsd(t2, getConstantD64(maxn).getMem(), x86::kCmpNLT);
			cc->subsd(t1, getConstantD64(magic0).getMem());

			// Combine the result.
			if (dst.getId() != src.getId())
				cc->movsd(dst, src);

			cc->andpd(dst, t2);
			cc->andnpd(t2, t1);
			cc->orpd(dst, t2);

			return;
		}

		// The `roundeven()` function can be used to implement efficiently the
		// remaining rounding functions. The easiest are `floor()` and `ceil()`.
		if (op == kOpRound || op == kOpFloor || op == kOpCeil) {
			X86Xmm t1 = cc->newXmmSd();
			X86Xmm t2 = cc->newXmmSd();
			X86Xmm t3 = cc->newXmmSd();

			cc->movsd(t2, src);
			cc->movsd(t3, src);

			if (dst.getId() != src.getId())
				cc->movsd(dst, src);

			switch (op) {
			case kOpRound:
				cc->addsd(t2, getConstantD64(magic0).getMem());
				cc->addsd(t3, getConstantD64(magic1).getMem());

				cc->movsd(t1, src);
				cc->subsd(t2, getConstantD64(magic0).getMem());
				cc->subsd(t3, getConstantD64(magic1).getMem());

				cc->cmpsd(t1, getConstantD64(maxn).getMem(), x86::kCmpNLT);
				cc->maxsd(t2, t3);

				cc->andpd(dst, t1);
				break;

			case kOpFloor:
				cc->addsd(t2, getConstantD64(magic0).getMem());
				cc->movsd(t1, src);

				cc->subsd(t2, getConstantD64(magic0).getMem());
				cc->cmpsd(t1, getConstantD64(maxn).getMem(), x86::kCmpNLT);

				cc->cmpsd(t3, t2, x86::kCmpLT);
				cc->andpd(t3, getConstantD64AsPD(1.0).getMem());

				cc->andpd(dst, t1);
				cc->subpd(t2, t3);
				break;

			case kOpCeil:
				cc->addsd(t2, getConstantD64(magic0).getMem());
				cc->movsd(t1, src);

				cc->subsd(t2, getConstantD64(magic0).getMem());
				cc->cmpsd(t1, getConstantD64(maxn).getMem(), x86::kCmpNLT);

				cc->cmpsd(t3, t2, x86::kCmpNLE);
				cc->andpd(t3, getConstantD64AsPD(1.0).getMem());

				cc->andpd(dst, t1);
				cc->addpd(t2, t3);
				break;
			}

			cc->andnpd(t1, t2);
			cc->orpd(dst, t1);

			return;
		}

		if (op == kOpTrunc) {
			X86Xmm t1 = cc->newXmmSd();
			X86Xmm t2 = cc->newXmmSd();
			X86Xmm t3 = cc->newXmmSd();

			cc->movsd(t2, getConstantU64(ASMJIT_UINT64_C(0x7FFFFFFFFFFFFFFF)).getMem());
			cc->andpd(t2, src);

			if (dst.getId() != src.getId())
				cc->movsd(dst, src);

			cc->movsd(t1, t2);
			cc->addsd(t2, getConstantD64(magic0).getMem());
			cc->movsd(t3, t1);

			cc->subsd(t2, getConstantD64(magic0).getMem());
			cc->cmpsd(t1, getConstantD64(maxn).getMem(), x86::kCmpNLT);

			cc->cmpsd(t3, t2, x86::kCmpLT);
			cc->orpd(t1, getConstantU64AsPD(ASMJIT_UINT64_C(0x8000000000000000)).getMem());
			cc->andpd(t3, getConstantD64AsPD(1.0).getMem());

			cc->andpd(dst, t1);
			cc->subpd(t2, t3);

			cc->andnpd(t1, t2);
			cc->orpd(dst, t1);

			return;
		}

		inlineCallAbstract(dst, &src, 2, op, isComplex, returnsComplex);
	}

	void JitCompiler::inlineCallAbstract(const X86Xmm& dst, const X86Xmm* args, uint32_t count, uint32_t op, bool paramsAreComplex, bool returnsComplex) {
		if (paramsAreComplex) {
			if (returnsComplex)
				inlineCallComplex(dst, args, count, JitUtils::getFuncByOp(op));
			else
				inlineCallCRetD(dst, args, count, JitUtils::getFuncByOp(op));
		}
		else {
			if (returnsComplex)
				inlineCallDRetC(dst, args, count, JitUtils::getFuncByOp(op));
			else
				inlineCall(dst, args, count, JitUtils::getFuncByOp(op));

		}
	}

	void JitCompiler::inlineCall(const X86Xmm& dst, const X86Xmm* args, uint32_t count, void* fn) {
		uint32_t i;

		// Use function builder to build a function prototype.
		FuncSignatureX signature;
		signature.setRetT<double>();

		for (i = 0; i < count; i++)
			signature.addArgT<double>();

		// Create the function call.
		CCFuncCall* ctx = cc->call((uint64_t)fn, signature);
		ctx->setRet(0, dst);

		for (i = 0; i < count; i++)
			ctx->setArg(static_cast<uint32_t>(i), args[i]);
	}

	void JitCompiler::inlineCallDRetC(const X86Xmm& dst, const X86Xmm* args, uint32_t count, void* fn) {
		size_t i;

		uint32_t length = count * sizeof(double);

		// copy the parameters to Memory.
		X86Mem stack = cc->newStack(length, sizeof(double));
		X86Gp dataPointerReg = cc->newUIntPtr();
		cc->lea(dataPointerReg, stack);

		// allocate Memory for the return
		X86Mem ret = cc->newStack(16, 16);
		X86Gp retReg = cc->newUIntPtr();
		cc->lea(retReg, ret);

		for (i = 0; i < count; i++)
		{
			cc->movsd(stack, args[i]);
			stack.addOffset(sizeof(double));
		}
		stack.resetOffset();


		X86Gp _ptr = cc->newUIntPtr();

		if (cc->getArchInfo().is64Bit()) {
			X86Mem function = cc->newUInt64Const(kAstScopeGlobal, (uint64_t)fn);
			cc->mov(_ptr, function);
		}
		else {
			X86Mem function = cc->newUInt32Const(kAstScopeGlobal, uint32_t((uint64_t)fn));
			cc->mov(_ptr, function);
		}


		// Use function builder to build a function prototype.
		FuncSignatureX signature;
		signature.setRetT<void>();
		signature.addArgT<TypeId::UIntPtr>(); // function-pointer
		signature.addArgT<TypeId::UIntPtr>(); // pointer to the return
		signature.addArgT<TypeId::UIntPtr>(); // parameters

		CCFuncCall* ctx;
		// Create the function call.
		ctx = cc->call((uint64_t)mpWrapDoubleC, signature);

		ctx->setArg(0, _ptr);
		ctx->setArg(1, retReg);
		ctx->setArg(2, dataPointerReg);

		cc->movapd(dst, ret);
	}

	//! Calls a function with complex arguments and non-complex returns.
	void JitCompiler::inlineCallCRetD(const X86Xmm& dst, const X86Xmm* args, const uint32_t count, void* fn) {
		size_t i;

		uint32_t length = count * 16;

		// copy the data to Memory.
		X86Mem stack = cc->newStack(length, 16);
		X86Gp ind = cc->newUIntPtr();
		cc->lea(ind, stack);

		for (i = 0; i < count; i++)
		{
			cc->movapd(stack, args[i]);
			stack.addOffset(16);
		}
		stack.resetOffset();

		X86Gp _ptr = cc->newUIntPtr();

		if (cc->getArchInfo().is64Bit()) {
			X86Mem function = cc->newUInt64Const(kAstScopeGlobal, (uint64_t)fn);
			cc->mov(_ptr, function);
		}
		else {
			X86Mem function = cc->newUInt32Const(kAstScopeGlobal, uint32_t((uint64_t)fn));
			cc->mov(_ptr, function);
		}


		// Use function builder to build a function prototype.
		FuncSignatureX signature;
		signature.setRetT<double>();
		signature.addArgT<TypeId::UIntPtr>(); // function-pointer
		signature.addArgT<TypeId::UIntPtr>(); // data

		CCFuncCall* ctx;
		// Create the function call.
		ctx = cc->call((uint64_t)mpWrapComplexD, signature);

		ctx->setRet(0, dst);
		ctx->setArg(0, _ptr);
		ctx->setArg(1, ind);

	}

	void JitCompiler::inlineCallComplex(const X86Xmm& dst, const X86Xmm* args, const uint32_t count, void* fn) {
		size_t i;

		uint32_t length = (count + 1) * 16;

		// copy the data to Memory.
		X86Mem stack = cc->newStack(length, 16);
		X86Gp ind = cc->newUIntPtr();
		cc->lea(ind, stack);

		for (i = 0; i < count; i++)
		{
			stack.addOffset(16);
			cc->movapd(stack, args[i]);
		}
		stack.resetOffset();

		X86Gp _ptr = cc->newUIntPtr();

		if (cc->getArchInfo().is64Bit()) {
			X86Mem function = cc->newUInt64Const(kAstScopeGlobal, (uint64_t)fn);
			cc->mov(_ptr, function);
		}
		else {
			X86Mem function = cc->newUInt32Const(kAstScopeGlobal, uint32_t((uint64_t)fn));
			cc->mov(_ptr, function);
		}


		// Use function builder to build a function prototype.
		FuncSignatureX signature;
		signature.setRetT<void>();
		signature.addArgT<TypeId::UIntPtr>();
		signature.addArgT<TypeId::UIntPtr>();

		CCFuncCall* ctx;
		// Create the function call.
		ctx = cc->call((uint64_t)mpWrapComplex, signature);

		ctx->setArg(0, _ptr);
		ctx->setArg(1, ind);

		cc->movapd(dst, stack);
	}

	void JitCompiler::prepareConstPool() {
		if (!constLabel.isValid()) {
			constLabel = cc->newLabel();

			CBNode* prev = cc->setCursor(functionBody);
			cc->lea(constPtr, x86::ptr(constLabel));
			if (prev != functionBody) cc->setCursor(prev);
		}
	}

	void JitCompiler::prepareConstPoolComp() {
		if (!constLabel.isValid()) {
			constLabel = cc->newLabel();

			CBNode* prev = cc->setCursor(functionBody);
			cc->lea(constPtr, x86::ptr(constLabel));
			if (prev != functionBody) cc->setCursor(prev);
		}
	}

	JitVar JitCompiler::getConstantU64(uint64_t value) {
		prepareConstPool();

		size_t offset;
		if (constPool.add(&value, sizeof(uint64_t), offset) != kErrorOk)
			return JitVar();

		return JitVar(x86::ptr(constPtr, static_cast<int>(offset)), JitVar::FLAG_NONE);
	}

	JitVar JitCompiler::getConstantU64Compl(uint64_t real, uint64_t imag) {
		prepareConstPoolComp();

		uint64_t value[2] = { real, imag };
		size_t offset;
		if (constPool.add(value, 2 * sizeof(uint64_t), offset) != kErrorOk)
			return JitVar();

		return JitVar(x86::ptr(constPtr, static_cast<int>(offset)), JitVar::FLAG_NONE);
	}

	JitVar JitCompiler::getConstantU64AsPD(uint64_t value) {
		prepareConstPool();

		size_t offset;
		Data128 vec = Data128::fromI64(value, 0);
		if (constPool.add(&vec, sizeof(Data128), offset) != kErrorOk)
			return JitVar();

		return JitVar(x86::ptr(constPtr, static_cast<int>(offset)), JitVar::FLAG_NONE);
	}

	JitVar JitCompiler::getConstantD64(double value) {
		DoubleBits bits;
		bits.d = value;
		return getConstantU64(bits.u);
	}

	JitVar JitCompiler::getConstantD64Compl(std::complex<double> value) {
		DoubleBitsComp bits;
		bits.d[0] = value.real();
		bits.d[1] = value.imag();
		return getConstantU64Compl(bits.u[0], bits.u[1]);
	}

	JitVar JitCompiler::getConstantD64AsPD(double value) {
		DoubleBits bits;
		bits.d = value;
		return getConstantU64AsPD(bits.u);
	}

	CompiledFunc mpCompileFunction(AstBuilder* ast, uint32_t options, OutputLog* log, bool b_complex) {
		StringLogger logger;

		CodeHolder code;
		code.init((jitGlobal.runtime.getCodeInfo()));

		X86Compiler c(&code);
		bool debugAsm = log != NULL && (options & kOptionDebugAsm) != 0;

		if (debugAsm) {
			logger.addOptions(Logger::kOptionBinaryForm | (b_complex ? Logger::kOptionImmExtended : 0));
			code.setLogger(&logger);
		}

		JitCompiler jitCompiler(ast->getHeap(), &c);
		if ((options & kOptionDisableSSE4_1) != 0)
			jitCompiler.enableSSE4_1 = false;

		jitCompiler.beginFunction();
		jitCompiler.compile(ast->getProgramNode(), ast->getRootScope(), ast->_numSlots, b_complex);
		jitCompiler.endFunction();

		c.finalize();

		CompiledFunc fn;
		jitGlobal.runtime.add(&fn, &code);

		if (debugAsm)
			log->log(OutputLog::kMessageAsm, 0, 0, logger.getString(), logger._stringBuilder.getLength());

		return fn;
	}

	void mpFreeFunction(void* fn) {
		jitGlobal.runtime.release(fn);
	}

} // mathpresso namespace
