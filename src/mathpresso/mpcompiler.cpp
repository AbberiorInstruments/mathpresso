// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include <mathpresso/mpast_p.h>
#include <mathpresso/mpcompiler_p.h>
#include <mathpresso/mpeval_p.h>
#include <mathpresso/mpoperation.h>


namespace mathpresso {
	using namespace asmjit;

	// ============================================================================
	// [mathpresso::JitGlobal]
	// ============================================================================

	struct JitGlobal {
		JitRuntime runtime;
	};
	static JitGlobal jitGlobal;


	void JitCompiler::beginFunction()
	{
		cc->addFunc(FuncSignature2<void, double*, double*>(CallConv::kIdHostCDecl));

		resultAddress = cc->newIntPtr("pResult");
		variablesAddress = cc->newIntPtr("pVariables");
		constPtr = cc->newIntPtr("pConst");

		cc->setArg(0, resultAddress);
		cc->setArg(1, variablesAddress);

		functionBody = cc->getCursor();
	}

	void JitCompiler::endFunction()
	{
		cc->endFunc();

		if (constLabel.isValid())
			cc->embedConstPool(constLabel, constPool);
	}

	JitVar JitCompiler::copyVar(const JitVar& other, uint32_t flags)
	{
		JitVar v(cc->newXmmSd(), flags);
		cc->emit(other.isXmm() ? X86Inst::kIdMovapd : X86Inst::kIdMovsd,
			v.getXmm(), other.getOperand());
		return v;
	}

	JitVar JitCompiler::writableVar(const JitVar& other)
	{
		if (other.isMem() || other.isRO())
			return copyVar(other, other.flags & ~JitVar::FLAG_RO);
		else
			return other;
	}

	// copies a non-complex var to a register, if necessary
	JitVar JitCompiler::registerVar(const JitVar& other)
	{
		if (other.isMem())
			return copyVar(other, other.flags);
		else
			return other;
	}

	JitVar JitCompiler::copyVarComplex(const JitVar& other, uint32_t flags)
	{
		JitVar v(cc->newXmmPd(), flags);

		cc->emit(X86Inst::kIdMovupd, v.getXmm(), other.getOperand());
		return v;
	}

	JitVar JitCompiler::writableVarComplex(const JitVar& other)
	{
		if (other.isMem() || other.isRO())
			return copyVarComplex(other, other.flags & ~JitVar::FLAG_RO);
		else
			return other;
	}

	// copies a var to a register, and makes sure its a complex value.
	JitVar JitCompiler::registerVarComplex(const JitVar& other, bool otherIsNonComplex)
	{
		if (otherIsNonComplex)
			return registerVarAsComplex(other);
		else if (other.isMem())
			return copyVarComplex(other, JitVar::FLAGS::FLAG_NONE);
		else
			return other;
	}

	// makes a non-complex var complex
	JitVar JitCompiler::registerVarAsComplex(const JitVar& other)
	{
		JitVar v(cc->newXmmPd(), other.flags);
		cc->pxor(v.getXmm(), v.getXmm());
		if (other.isMem())
		{
			cc->movlpd(v.getXmm(), other.getMem());
		}
		else
		{
			cc->movsd(v.getXmm(), other.getXmm());
		}
		return v;
	}

	//! Compiles an AstBlock into assembler.
	//! NOTE: use beginFunction() before and endFunction() after calling this.
	void JitCompiler::compile(AstBlock* node, AstScope* rootScope, uint32_t numSlots, bool b_complex)
	{
		// Create Definitions for the Variables and add them as JitVar
		if (numSlots != 0)
		{
			varSlots = static_cast<JitVar*>(heap->alloc(sizeof(JitVar) * numSlots));
			if (varSlots == nullptr) return;

			for (uint32_t i = 0; i < numSlots; i++)
				varSlots[i] = JitVar();
		}

		// Result of the function or NaN. Here the AST is compiled.
		JitVar result = onBlock(node);

		// Write altered global variables.
		{
			AstSymbolHashIterator it(rootScope->getSymbols());
			while (it.has())
			{
				AstSymbol* sym = it.get();
				if (sym->isGlobal() && sym->isAltered())
				{
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
				var = registerVarComplex(getConstantD64(mpGetNan())).getXmm();
			else
				var = registerVarComplex(result, !node->returnsComplex()).getXmm();
			cc->movupd(x86::ptr(resultAddress), var);
		}

		// Release the Space allocated for the variables
		if (numSlots != 0)
			heap->release(varSlots, sizeof(JitVar) * numSlots);
	}

	JitVar JitCompiler::onNode(AstNode* node)
	{
		switch (node->getNodeType())
		{
		case kAstNodeBlock: return onBlock(static_cast<AstBlock*>(node));
		case kAstNodeVar: return onVar(static_cast<AstVar*>(node));
		case kAstNodeImm: return onImm(static_cast<AstImm*>(node));
		case kAstNodeVarDecl:
		case kAstNodeUnaryOp: 
		case kAstNodeBinaryOp:
		case kAstNodeTernaryOp:
		case kAstNodeCall: 
			return node->_mpOp->compile(this, node);

		default:
			MATHPRESSO_ASSERT_NOT_REACHED();
			return JitVar();
		}
	}

	JitVar JitCompiler::onBlock(AstBlock* node)
	{
		JitVar result;
		size_t i, len = node->getLength();

		for (i = 0; i < len; i++)
			result = onNode(node->getAt(i));

		// Return the last result (or no result if the block is empty).
		return result;
	}

	JitVar JitCompiler::onVar(AstVar* node)
	{
		AstSymbol* sym = node->getSymbol();
		uint32_t slotId = sym->getVarSlotId();
		bool b_complex = node->returnsComplex();

		JitVar result = varSlots[slotId];
		if (result.isNone())
		{
			if (sym->isGlobal())
			{
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
				result = getConstantD64(mpGetNan());
				varSlots[slotId] = result;
			}
		}

		return result;
	}

	JitVar JitCompiler::onImm(AstImm* node)
	{
		if (node->returnsComplex())
			return getConstantD64(node->getValueCplx());
		else
			return getConstantD64(node->getValue());
	}

	void JitCompiler::inlineCallDRetD(const X86Xmm& dst, const X86Xmm* args, size_t count, void* fn)
	{

#ifdef _REALREWORK
		// Use function builder to build a function prototype.
		X86Mem stack(cc->newStack(count * sizeof(double), sizeof(double)));
		X86Gp dataPointerReg(cc->newUIntPtr());
		cc->lea(dataPointerReg, stack);

		for (size_t i = 0; i < count; i++)
		{
			cc->movsd(stack, args[i]);
			stack.addOffset(sizeof(double));
		}

		FuncSignatureX signature;
		signature.setRetT<double>();
		signature.addArgT<TypeId::UIntPtr>(); // parameters

		// Create the function call.
		CCFuncCall* ctx = cc->call(reinterpret_cast<uint64_t>(fn), signature);
		ctx->setRet(0, dst);

		ctx->setArg(0, dataPointerReg);

#else

		// Use function builder to build a function prototype.
		FuncSignatureX signature;
		signature.setRetT<double>();

		for (size_t i = 0; i < count; i++)
			signature.addArgT<double>();

		//Create the function call.
		CCFuncCall* ctx = cc->call(reinterpret_cast<uint64_t>(fn), signature);
		ctx->setRet(0, dst);

		for (size_t i = 0; i < count; i++)
			ctx->setArg(static_cast<uint32_t>(i), args[i]);

#endif // _REALREWORK
	}

	void JitCompiler::inlineCallDRetC(const X86Xmm& dst, const X86Xmm* args, size_t count, void* fn)
	{
		// copy the parameters to Memory.
		X86Mem stack(cc->newStack(uint32_t(count * sizeof(double)), sizeof(double)));
		X86Gp dataPointerReg(cc->newUIntPtr());
		cc->lea(dataPointerReg, stack);

		// allocate Memory for the return
		X86Mem ret(cc->newStack(16, 16));
		X86Gp retReg(cc->newUIntPtr());
		cc->lea(retReg, ret);

		for (size_t i = 0; i < count; i++)
		{
			cc->movsd(stack, args[i]);
			stack.addOffset(sizeof(double));
		}

		// Use function builder to build a function prototype.
		FuncSignatureX signature;
		signature.setRetT<void>();
		signature.addArgT<TypeId::UIntPtr>(); // function-pointer
		signature.addArgT<TypeId::UIntPtr>(); // pointer to the return
		signature.addArgT<TypeId::UIntPtr>(); // parameters

		// Create the function call.
		CCFuncCall* ctx = cc->call(reinterpret_cast<uint64_t>(mpWrapDtoC), signature);

		ctx->setArg(0, imm_u(reinterpret_cast<uint64_t>(fn)));
		ctx->setArg(1, dataPointerReg);
		ctx->setArg(2, retReg);

		cc->movapd(dst, ret);
	}

	//! Calls a function with complex arguments and non-complex returns.
	void JitCompiler::inlineCallCRetD(const X86Xmm& dst, const X86Xmm* args, const size_t count, void* fn)
	{
		// copy the data to Memory.
		X86Mem stack(cc->newStack(uint32_t(count) * 16, 16));
		X86Gp dataPointerReg(cc->newUIntPtr());
		cc->lea(dataPointerReg, stack);

		for (size_t i = 0; i < count; i++)
		{
			cc->movapd(stack, args[i]);
			stack.addOffset(16);
		}
		stack.resetOffset();

		// Use function builder to build a function prototype.
		FuncSignatureX signature;
		signature.setRetT<double>();
		signature.addArgT<TypeId::UIntPtr>(); // data

		CCFuncCall* ctx = cc->call(reinterpret_cast<uint64_t>(fn), signature);

		ctx->setRet(0, dst);
		ctx->setArg(0, dataPointerReg);

	}

	void JitCompiler::inlineCallCRetC(const X86Xmm& dst, const X86Xmm* args, size_t count, void* fn)
	{
		// copy the data to Memory.
		X86Mem stack(cc->newStack((uint32_t(count) + 1) * 16, 16));
		X86Gp dataPointerReg(cc->newUIntPtr());
		cc->lea(dataPointerReg, stack);

		for (size_t i = 0; i < count; i++)
		{
			stack.addOffset(16);
			cc->movapd(stack, args[i]);
		}
		stack.resetOffset();

		// Use function builder to build a function prototype.
		FuncSignatureX signature;
		signature.setRetT<void>();
		signature.addArgT<TypeId::UIntPtr>();
		signature.addArgT<TypeId::UIntPtr>();

		// Create the function call.
		CCFuncCall* ctx = cc->call(reinterpret_cast<uint64_t>(mpWrapCtoC), signature);

		ctx->setArg(0, imm_u(reinterpret_cast<uint64_t>(fn)));
		ctx->setArg(1, dataPointerReg);

		cc->movapd(dst, stack);
	}

	void JitCompiler::prepareConstPool()
	{
		if (!constLabel.isValid())
		{
			constLabel = cc->newLabel();

			CBNode* prev = cc->setCursor(functionBody);
			cc->lea(constPtr, x86::ptr(constLabel));
			if (prev != functionBody) cc->setCursor(prev);
		}
	}

	JitVar JitCompiler::getConstantU64(uint64_t value)
	{
		prepareConstPool();

		size_t offset;
		if (constPool.add(&value, sizeof(uint64_t), offset) != kErrorOk)
			return JitVar();

		return JitVar(x86::ptr(constPtr, static_cast<int>(offset)), JitVar::FLAG_NONE);
	}

	JitVar JitCompiler::getConstantU64(uint64_t real, uint64_t imag)
	{
		prepareConstPool();

		uint64_t value[2] = { real, imag };
		size_t offset;
		if (constPool.add(value, 2 * sizeof(uint64_t), offset) != kErrorOk)
			return JitVar();

		return JitVar(x86::ptr(constPtr, static_cast<int>(offset)), JitVar::FLAG_NONE);
	}

	JitVar JitCompiler::getConstantU64AsPD(uint64_t value)
	{
		prepareConstPool();

		size_t offset;
		Data128 vec = Data128::fromI64(value, 0);
		if (constPool.add(&vec, sizeof(Data128), offset) != kErrorOk)
			return JitVar();

		return JitVar(x86::ptr(constPtr, static_cast<int>(offset)), JitVar::FLAG_NONE);
	}

	JitVar JitCompiler::getConstantD64(double value)
	{
		DoubleBits bits;
		bits.d = value;
		return getConstantU64(bits.u);
	}

	JitVar JitCompiler::getConstantD64(std::complex<double> value)
	{
		DoubleBitsComp bits;
		bits.d[0] = value.real();
		bits.d[1] = value.imag();
		return getConstantU64(bits.u[0], bits.u[1]);
	}

	JitVar JitCompiler::getConstantD64AsPD(double value)
	{
		DoubleBits bits;
		bits.d = value;
		return getConstantU64AsPD(bits.u);
	}

	CompiledFunc mpCompileFunction(AstBuilder* ast, uint32_t options, OutputLog* log, const Operations * ops, bool b_complex)
	{
		StringLogger logger;

		CodeHolder code;
		code.init((jitGlobal.runtime.getCodeInfo()));

		X86Compiler c(&code);
		bool debugAsm = log != nullptr && (options & kOptionDebugAsm) != 0;

		if (debugAsm)
		{
			logger.addOptions(Logger::kOptionBinaryForm | (b_complex ? Logger::kOptionImmExtended : 0));
			code.setLogger(&logger);
		}

		JitCompiler jitCompiler(ast->getHeap(), &c, ops);
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

	void mpFreeFunction(void* fn)
	{
		jitGlobal.runtime.release(fn);
	}

} // mathpresso namespace
