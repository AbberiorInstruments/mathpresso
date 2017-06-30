// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

#include "mpoperation_p.h"
#include "asmjit\x86\x86operand.h"


namespace mathpresso {


	JitVar OperationGeneric::compileToASM(JitCompiler* jc, AstNode * node) {
		// TODO: error handling?
		// TODO: switch to other fnXtoX, if initial not found? not necessary, because of optimize?
		// TODO: compile directly to assembler if asmXToX is found?

		asmjit::X86Xmm result = node->returnsComplex() ? jc->cc->newXmmPd() : jc->cc->newXmmSd();
		asmjit::X86Xmm args[8];

		if (!node->takesComplex())
		{
			for (size_t i = 0; i < numArgs_; i++)
			{
				args[i] = jc->registerVar(jc->onNode(node->getAt(i))).getXmm();
			}

			if (!node->returnsComplex())
			{
				if (fnDToD)
				{
					jc->inlineCall(result, args, numArgs_, fnDToD);
				}
				else
				{
					// do some throwing?
					return JitVar();
				}
			}
			else
			{
				if (fnDToC)
				{
					jc->inlineCallDRetC(result, args, numArgs_, fnDToC);
				}
				else
				{
					// do some throwing?
					return JitVar();
				}
			}
		}
		else
		{
			for (size_t i = 0; i < numArgs_; i++)
			{
				args[i] = jc->registerVarComplex(jc->onNode(node->getAt(i)), !node->getAt(i)->returnsComplex()).getXmm();
			}

			if (!node->takesComplex())
			{
				if (fnCToD)
				{
					jc->inlineCallCRetD(result, args, numArgs_, fnCToD);
				}
				else
				{
					// do some throwing?
					return JitVar();
				}
			}
			else
			{
				if (fnCToC)
				{
					jc->inlineCallComplex(result, args, numArgs_, fnCToC);
				}
				else
				{
					// do some throwing?
					return JitVar();
				}
			}
		}

		return JitVar(result, JitVar::FLAG_NONE);
	}

	Error OperationGeneric::optimize(AstOptimizer * opt, AstNode * node) {

		uint32_t count = node->getLength();

		bool b_need_cplx = false, b_all_imm = true;
		for (size_t i = 0; i < count; i++)
		{
			opt->onNode(node->getAt(i));
			b_need_cplx |= node->getAt(i)->returnsComplex();
			b_all_imm |= node->getAt(i)->isImm();
		}

		// set flags according to the available functions.
		if (b_need_cplx)
		{
			node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
			if (fnCToD)
			{
				node->removeNodeFlags(AstNodeFlags::kAstReturnsComplex);
			}
			else if (fnCToC)
			{
				node->addNodeFlags(AstNodeFlags::kAstReturnsComplex);
			}
			else
			{
				return ErrorCode::kErrorInvalidState;
			}
		}
		else
		{
			if (fnDToD)
			{
				node->removeNodeFlags(AstNodeFlags::kAstReturnsComplex | AstNodeFlags::kAstTakesComplex);
			}
			else if (fnCToD)
			{
				node->removeNodeFlags(AstNodeFlags::kAstTakesComplex);
				node->addNodeFlags(AstNodeFlags::kAstReturnsComplex);
			}
			else if (fnDToC)
			{
				node->removeNodeFlags(AstNodeFlags::kAstReturnsComplex);
				node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
			}
			else if (fnCToC)
			{
				node->addNodeFlags(AstNodeFlags::kAstReturnsComplex | AstNodeFlags::kAstTakesComplex);
			}
			else
			{
				return ErrorCode::kErrorInvalidState;
			}
		}

		if (b_all_imm) {
			AstImm* ret = static_cast<AstImm*>(node->getAt(0));
			if (node->takesComplex())
			{
				std::complex<double> args[8];
				for (size_t i = 0; i < numArgs_; i++)
				{
					args[i] = (static_cast<AstImm*>(node->getAt(i)))->getValueCplx();
				}

				if (node->returnsComplex())
				{
					ret->setValue(evaluateCtoC(args));
				}
				else
				{
					ret->setValue(evaluateCtoD(args));
				}
			}
			else
			{
				double args[8];
				for (size_t i = 0; i < numArgs_; i++)
				{
					args[i] = (static_cast<AstImm*>(node->getAt(i)))->getValue();
				}

				if (node->returnsComplex())
				{
					ret->setValue(evaluateDtoC(args));
				}
				else
				{
					ret->setValue(evaluateDtoD(args));
				}
			}
		}

		return ErrorCode::kErrorOk;
	}


	std::complex<double> OperationGeneric::evaluateCtoC(std::complex<double>* args) {
		if (fnDToD)
		{
			return ((mpFuncpCtoC)fnCToC)(args);
		}
		else
		{
			// throw?
			//return ->_errorReporter->onError(kErrorInvalidState, node->getPosition(),
			//"Invalid binary operation '%s'.", name);
			return std::complex<double>(0, 0);
		}
	}

	std::complex<double> OperationGeneric::evaluateDtoC(double * args) {
		if (fnDToC)
		{
			return ((mpFuncpDtoC)fnDToC)(args);
		}
		else
		{
			// throw?
			//return ->_errorReporter->onError(kErrorInvalidState, node->getPosition(),
			//"Invalid binary operation '%s'.", name);
			return std::complex<double>(0, 0);
		}
	}

	double OperationGeneric::evaluateCtoD(std::complex<double>* args) {
		if (fnCToD)
		{
			return ((mpFuncpCtoD)fnCToD)(args);
		}
		else
		{
			// throw?
			//return ->_errorReporter->onError(kErrorInvalidState, node->getPosition(),
			//"Invalid binary operation '%s'.", name);
			return 0;
		}
	}

	double OperationGeneric::evaluateDtoD(double * args) {
		if (fnDToD)
		{
			switch (numArgs_)
			{
			case 0: return ((Arg0Func)fnDToD)();
			case 1: return ((Arg1Func)fnDToD)(args[0]);
			case 2: return ((Arg2Func)fnDToD)(args[0], args[1]);
			case 3: return ((Arg3Func)fnDToD)(args[0], args[1], args[2]);
			case 4: return ((Arg4Func)fnDToD)(args[0], args[1], args[2], args[3]);
			case 5: return ((Arg5Func)fnDToD)(args[0], args[1], args[2], args[3], args[4]);
			case 6: return ((Arg6Func)fnDToD)(args[0], args[1], args[2], args[3], args[4], args[5]);
			case 7: return ((Arg7Func)fnDToD)(args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
			case 8: return ((Arg8Func)fnDToD)(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]);
			default:
				// throw something?
				return 0;
			}
		}
		else
		{
			return 0;
			// throw something?
		}
	}
}
