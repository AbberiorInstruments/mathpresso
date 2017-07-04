// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

#include "mpoperation_p.h"
#include "asmjit\x86\x86operand.h"


namespace mathpresso {

	JitVar MpOperationFunc::compile(JitCompiler* jc, AstNode * node) 
	{
		asmjit::X86Xmm result = node->returnsComplex() ? jc->cc->newXmmPd() : jc->cc->newXmmSd();
		asmjit::X86Xmm args[8];

		if (!node->takesComplex())
		{
			for (size_t i = 0; i < nargs_; i++)
			{
				args[i] = jc->registerVar(jc->onNode(node->getAt(i))).getXmm();
			}

			bool returnsComplex = (flags_ & OpFlagDReturnsC) != 0;
			if (node->returnsComplex() != returnsComplex  || !fnD_)
			{
				// Should never happen, as the optimizer should have taken care of that. Remove later
				throw std::runtime_error("Implementation error!");
			}

			if (returnsComplex)
			{
				jc->inlineCallDRetC(result, args, nargs_, fnD_);
			}
			else
			{
				jc->inlineCallDRetD(result, args, nargs_, fnD_);
			}
		}
		else
		{
			for (size_t i = 0; i < nargs_; i++)
			{
				args[i] = jc->registerVarComplex(jc->onNode(node->getAt(i)), !node->getAt(i)->returnsComplex()).getXmm();
			}

			bool returnsComplex = (flags_ & OpFlagCReturnsD) == 0;
			if (node->returnsComplex() != returnsComplex || !fnC_)
			{
				// Should never happen, as the optimizer should have taken care of that. Remove later
				throw std::runtime_error("Implementation error!");
			}

			if (returnsComplex)
			{
				jc->inlineCallCRetC(result, args, nargs_, fnC_);
			}
			else
			{
				jc->inlineCallCRetD(result, args, nargs_, fnC_);
			}
		}

		return JitVar(result, JitVar::FLAG_NONE);
	}

	Error MpOperationFunc::optimize(AstOptimizer * opt, AstNode * node) 
	{
		size_t count = node->getLength();
		bool b_need_cplx = false;
		bool b_all_imm = true;

		for (size_t i = 0; i < count; i++)
		{
			opt->onNode(node->getAt(i));
			b_need_cplx |= node->getAt(i)->returnsComplex();
			b_all_imm |= node->getAt(i)->isImm();
		}

		bool b_returns_complex;

		// set flags according to the available functions.
		if (b_need_cplx)
		{
			if (!fnC_)
				return ErrorCode::kErrorInvalidArgument;

			node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
			b_returns_complex = flags_ & OpFlagCReturnsD;
		}
		else
		{
			if (fnD_)
			{
				node->removeNodeFlags(AstNodeFlags::kAstTakesComplex);
				b_returns_complex = flags_ & OpFlagDReturnsC;
			}
			else if (fnC_)
			{
				node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
				b_returns_complex = (flags_ & OpFlagCReturnsD) == 0;
			}
			else
			{
				return ErrorCode::kErrorSymbolNotFound;
			}
		}
					
		if (b_returns_complex)
		{
			node->addNodeFlags(AstNodeFlags::kAstReturnsComplex);
		}
		else
		{
			node->removeNodeFlags(AstNodeFlags::kAstReturnsComplex);
		}


		if (b_all_imm) 
		{
			AstImm* ret = static_cast<AstImm*>(node->getAt(0));
			if (node->takesComplex())
			{
				std::complex<double> args[8];
				for (size_t i = 0; i < nargs_; i++)
				{
					args[i] = (static_cast<AstImm*>(node->getAt(i)))->getValueCplx();
				}

				if (b_returns_complex)
				{
					ret->setValue(((mpFuncpCtoC)fnC_)(args));
				}
				else
				{
					ret->setValue(((mpFuncpCtoD)fnC_)(args));
				}
			}
			else
			{
				double args[8];
				for (size_t i = 0; i < nargs_; i++)
				{
					args[i] = (static_cast<AstImm*>(node->getAt(i)))->getValue();
				}

				if (b_returns_complex)
				{
					ret->setValue(((mpFuncpDtoC)fnD_)(args));
				}
				else
				{
					ret->setValue(evaluateDRetD(args));
				}
			}
		}

		return ErrorCode::kErrorOk;
	}


	
	double MpOperationFunc::evaluateDRetD(double * args) 
	{
		switch (nargs_)
		{
		case 0: return ((Arg0Func)fnD_)();
		case 1: return ((Arg1Func)fnD_)(args[0]);
		case 2: return ((Arg2Func)fnD_)(args[0], args[1]);
		case 3: return ((Arg3Func)fnD_)(args[0], args[1], args[2]);
		case 4: return ((Arg4Func)fnD_)(args[0], args[1], args[2], args[3]);
		case 5: return ((Arg5Func)fnD_)(args[0], args[1], args[2], args[3], args[4]);
		case 6: return ((Arg6Func)fnD_)(args[0], args[1], args[2], args[3], args[4], args[5]);
		case 7: return ((Arg7Func)fnD_)(args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
		case 8: return ((Arg8Func)fnD_)(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]);
		default:
			throw std::runtime_error("Too many arguments.");
		}
	}
}
