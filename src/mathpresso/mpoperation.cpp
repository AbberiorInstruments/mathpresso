// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

#include "mpoperation_internal_p.h"
#include "mpast_p.h"
#include "mpcompiler_p.h"
#include "mpoptimizer_p.h"
#include "asmjit/x86/x86operand.h"
#include "asmjit/x86/x86inst.h"

#include <complex>

namespace mathpresso {

	uint32_t addBuiltinMpObjects(Context * ctx)
	{
		ctx->addObject("+", std::make_shared<MpOperationAdd>());
		ctx->addObject("-", std::make_shared<MpOperationSub>());
		ctx->addObject("*", std::make_shared<MpOperationMul>());
		ctx->addObject("/", std::make_shared<MpOperationDiv>());
		ctx->addObject("==", std::make_shared<MpOperationEq>());
		ctx->addObject("!=", std::make_shared<MpOperationNe>());
		ctx->addObject(">=", std::make_shared<MpOperationGe>());
		ctx->addObject(">", std::make_shared<MpOperationGt>());
		ctx->addObject("<=", std::make_shared<MpOperationLe>());
		ctx->addObject("<", std::make_shared<MpOperationLt>());
		ctx->addObject("?", std::make_shared<MpOperationTernary>());
		ctx->addObject("=", std::make_shared<MpOperationAssignment>());
		ctx->addObject("isfinite", std::make_shared<MpOperationIsFinite>());
		ctx->addObject("isinf", std::make_shared<MpOperationIsInfinite>());
		ctx->addObject("isnan", std::make_shared<MpOperationIsNan>());
		ctx->addObject("real", std::make_shared<MpOperationGetReal>());
		ctx->addObject("imag", std::make_shared<MpOperationGetImag>());
		ctx->addObject("min", std::make_shared<MpOperationMin>());
		ctx->addObject("max", std::make_shared<MpOperationMax>());
		ctx->addObject("=", std::make_shared<MpOperationAssignment>());
		ctx->addObject("sin", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::sin));
		ctx->addObject("cos", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::cos));
		ctx->addObject("tan", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::tan));
		ctx->addObject("sinh", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::sinh));
		ctx->addObject("cosh", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::cosh));
		ctx->addObject("tanh", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::tanh));
		ctx->addObject("asin", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::asin));
		ctx->addObject("acos", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::acos));
		ctx->addObject("atan", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::atan));
		ctx->addObject("%", std::make_shared<MpOperationModulo>());
		ctx->addObject("-", std::make_shared<MpOperationNeg>());
		ctx->addObject("!", std::make_shared<MpOperationNot>());
		ctx->addObject("sqrt", std::make_shared<MpOperationSqrt>());
		ctx->addObject("sqrtc", std::make_shared<MpOperationSqrtC>());
		ctx->addObject("conjug", std::make_shared<MpOperationConjug>());
		ctx->addObject("avg", std::make_shared<MpOperationAvg>());
		ctx->addObject("abs", std::make_shared<MpOperationAbs>());
		ctx->addObject("recip", std::make_shared<MpOperationRecip>());
		ctx->addObject("signbit", std::make_shared<MpOperationSignBit>());
		ctx->addObject("copysign", std::make_shared<MpOperationCopySign>());
		ctx->addObject("round", std::make_shared<MpOperationRound>());
		ctx->addObject("roundeven", std::make_shared<MpOperationRoundEven>());
		ctx->addObject("floor", std::make_shared<MpOperationFloor>());
		ctx->addObject("ceil", std::make_shared<MpOperationcCeil>());
		ctx->addObject("frac", std::make_shared<MpOperationFrac>());
		ctx->addObject("trunc", std::make_shared<MpOperationTrunc>());
		ctx->addObject("log", std::make_shared<MpOperationLog>());
		ctx->addObject("log2", std::make_shared<MpOperationLog2>());
		ctx->addObject("log10", std::make_shared<MpOperationLog10>());
		ctx->addObject("exp", std::make_shared<MpOperationExp>());
		ctx->addObject("pow", std::make_shared<MpOperationPow>());
		ctx->addObject("atan2", std::make_shared<MpOperationAtan2>());
		ctx->addObject("hypot", std::make_shared<MpOperationHypot>());
		ctx->addObject("_none_", std::make_shared<MpOperationFunc>(0, MpOperationFlags::OpFlagNone, nullptr, nullptr));
		return 0;
	}


	// MpOperationFunc
	JitVar MpOperationFunc::compile(JitCompiler* jc, AstNode * node)
	{
		asmjit::X86Xmm result = node->returnsComplex() ? jc->cc->newXmmPd() : jc->cc->newXmmSd();
		asmjit::X86Xmm args[8];

		if (!node->takesComplex())
		{
			bool returnsComplex = hasFlag(MpOperationFlags::OpFlagDReturnsC);
			if (node->returnsComplex() != returnsComplex || !fnD_)
			{
				// Should never happen, as the optimizer should have taken care of that. Remove later
				throw std::runtime_error("Implementation error!");
			}

			for (size_t i = 0; i < nargs_; i++)
			{
				args[i] = jc->registerVar(jc->onNode(node->getAt(i))).getXmm();
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
			bool returnsComplex = !hasFlag(MpOperationFlags::OpFlagCReturnsD);
			if (node->returnsComplex() != returnsComplex || !fnC_)
			{
				// Should never happen, as the optimizer should have taken care of that. Remove later
				throw std::runtime_error("Implementation error!");
			}

			for (size_t i = 0; i < nargs_; i++)
			{
				args[i] = jc->registerVarComplex(jc->onNode(node->getAt(i)), !node->getAt(i)->returnsComplex()).getXmm();
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

		// Gather Information about the child-nodes.
		for (size_t i = 0; i < count; i++)
		{
			MATHPRESSO_PROPAGATE(opt->onNode(node->getAt(i)));
			b_need_cplx |= node->getAt(i)->returnsComplex();
			b_all_imm &= node->getAt(i)->isImm();
		}

		bool b_returns_complex;

		// set flags according to the available functions.
		if (b_need_cplx)
		{
			if (hasFlag(MpOperationFlags::OpHasNoComplex))
				return opt->_errorReporter->onError(kErrorInvalidArgument, node->getPosition(),
					"No complex function available.");

			node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
			b_returns_complex = !hasFlag(MpOperationFlags::OpFlagCReturnsD);
		}
		else
		{
			if (!hasFlag(MpOperationFlags::OpHasNoReal))
			{
				node->removeNodeFlags(AstNodeFlags::kAstTakesComplex);
				b_returns_complex = hasFlag(MpOperationFlags::OpFlagDReturnsC);
			}
			else if (!hasFlag(MpOperationFlags::OpHasNoComplex))
			{
				node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
				b_returns_complex = !hasFlag(MpOperationFlags::OpFlagCReturnsD);
			}
			else
			{
				return opt->_errorReporter->onError(kErrorSymbolNotFound, node->getPosition(),
					"No appropriate function available");
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

		// optimize all-immediate calls:
		if (b_all_imm)
		{
			AstImm* ret = opt->getAst()->newNode<AstImm>(0);
			if (node->takesComplex())
			{
				std::vector<std::complex<double>> args;
				for (size_t i = 0; i < nargs_; i++)
				{
					args.push_back((static_cast<AstImm*>(node->getAt(i)))->getValueCplx());
				}

				if (node->returnsComplex())
				{
					ret->setValue(evaluateCRetC(args.data()));
				}
				else
				{
					ret->setValue(evaluateCRetD(args.data()));
				}
			}
			else
			{
				std::vector<double> args;
				for (size_t i = 0; i < nargs_; i++)
				{
					args.push_back((static_cast<AstImm*>(node->getAt(i)))->getValue());
				}

				if (node->returnsComplex())
				{
					ret->setValue(evaluateDRetC(args.data()));
				}
				else
				{
					ret->setValue(evaluateDRetD(args.data()));
				}
			}
			node->getParent()->replaceNode(node, ret);
			opt->onNode(ret);

			opt->getAst()->deleteNode(node);
		}
		return optimizeSpecial(opt, node);
	}

	void MpOperationFunc::setFn(void* fn, bool isComplex)
	{
		if (isComplex)
		{
			fnC_ = fn;
			flags_ &= ~OpHasNoComplex;
		}
		else
		{
			fnD_ = fn;
			flags_ &= ~OpHasNoReal;
		}
	}

	double MpOperationFunc::evaluateDRetD(double * args)
	{
		if (!fnD_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncDtoD)fnD_)(args);

		/*switch (nargs_)
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
		}*/
	}

	std::complex<double> MpOperationFunc::evaluateDRetC(double * args)
	{
		if (!fnD_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpDtoC)fnD_)(args);
	}

	double MpOperationFunc::evaluateCRetD(std::complex<double>* args)
	{
		if (!fnC_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpCtoD)fnC_)(args);
	}

	std::complex<double> MpOperationFunc::evaluateCRetC(std::complex<double>* args)
	{
		if (!fnC_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpCtoC)fnC_)(args);
	}

	uint32_t MpOperationFunc::optimizeSpecial(AstOptimizer * opt, AstNode * node) { return mathpresso::ErrorCode::kErrorOk; }

	// MpOperationFuncAsm
	JitVar mathpresso::MpOperationFuncAsm::compile(JitCompiler * jc, AstNode * node)
	{
		if (!hasFlag(MpOperationFlags::OpFlagHasAsm))
		{
			return MpOperationFunc::compile(jc, node);
		}
		else
		{
			if (!node->takesComplex())
			{
				bool returnsComplex = hasFlag(MpOperationFlags::OpFlagDReturnsC);
				if (node->returnsComplex() != returnsComplex || !fnD_)
				{
					// Should never happen, as the optimizer should have taken care of that. Remove later
					throw std::runtime_error("Implementation error!");
				}

				std::vector<JitVar> args;
				for (size_t i = 0; i < node->getLength(); i++)
				{
					args.push_back(jc->registerVar(jc->onNode(node->getAt(i))));
				}

				return asmD_(jc, args.data());
			}
			else
			{
				bool returnsComplex = !hasFlag(MpOperationFlags::OpFlagCReturnsD);
				if (node->returnsComplex() != returnsComplex || !fnC_)
				{
					// Should never happen, as the optimizer should have taken care of that. Remove later
					throw std::runtime_error("Implementation error!");
				}

				std::vector<JitVar> args;
				for (size_t i = 0; i < node->getLength(); i++)
				{
					args.push_back(jc->registerVarComplex(jc->onNode(node->getAt(i))));
				}

				return asmC_(jc, args.data());
			}
		}
	}

	void MpOperationFuncAsm::setFnAsm(mpAsmFunc fn, bool isComplex)
	{
		if (isComplex)
		{
			asmC_ = fn;
		}
		else
		{
			asmD_ = fn;
		}
	}

	// MpOperationIsFinite
	JitVar MpOperationIsFinite::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var = jc->onNode(node->getAt(0));

		if (node->takesComplex())
		{
			var = jc->writableVarComplex(var);
			jc->cc->orpd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0x8000000000000000), MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
			jc->cc->cmppd(var.getXmm(), jc->getConstantD64(std::complex<double>(0.0, 0.0)).getMem(), int(asmjit::x86::kCmpLE));
			jc->cc->andpd(var.getXmm(), jc->getConstantD64(std::complex<double>(1.0, 1.0)).getMem());
		}
		else
		{
			var = jc->writableVar(var);
			jc->cc->orpd(var.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
			jc->cc->cmpsd(var.getXmm(), jc->getConstantU64(0).getMem(), int(asmjit::x86::kCmpLE));
			jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		}
		return var;
	}

	double MpOperationIsFinite::evaluateDRetD(double * args) { return std::isfinite(args[0]) ? 1.0 : 0.0; }

	std::complex<double> MpOperationIsFinite::evaluateCRetC(std::complex<double>* args)
	{
		return std::complex<double>(std::isfinite(args[0].real()) ? 1.0 : 0.0, std::isfinite(args[0].imag()) ? 1.0 : 0.0);
	}

	// MpOperationIsInFinite
	JitVar MpOperationIsInfinite::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var = jc->onNode(node->getAt(0));

		if (node->takesComplex())
		{
			var = jc->writableVarComplex(var);
			jc->cc->orpd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0x8000000000000000), MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
			jc->cc->cmppd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0xFFF0000000000000), MATHPRESSO_UINT64_C(0xFFF0000000000000)).getMem(), int(asmjit::x86::kCmpEQ));
			jc->cc->andpd(var.getXmm(), jc->getConstantD64(std::complex<double>(1.0, 1.0)).getMem());
		}
		else
		{
			var = jc->writableVar(var);
			jc->cc->orpd(var.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
			jc->cc->cmpsd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0xFFF0000000000000)).getMem(), int(asmjit::x86::kCmpEQ));
			jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		}
		return var;
	}

	double MpOperationIsInfinite::evaluateDRetD(double * args) { return std::isinf(args[0]) ? 1.0 : 0.0; }

	std::complex<double> MpOperationIsInfinite::evaluateCRetC(std::complex<double>* args)
	{
		return std::complex<double>(std::isinf(args[0].real()) ? 1.0 : 0.0, std::isinf(args[0].imag()) ? 1.0 : 0.0);
	}

	// MpOperationIsNan
	JitVar MpOperationIsNan::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var = jc->onNode(node->getAt(0));

		if (node->takesComplex())
		{
			var = jc->writableVarComplex(var);
			jc->cc->cmppd(var.getXmm(), var.getXmm(), int(asmjit::x86::kCmpEQ));
			jc->cc->andnpd(var.getXmm(), jc->getConstantD64(std::complex<double>(1.0, 1.0)).getMem());
		}
		else
		{
			var = jc->writableVar(var);
			jc->cc->cmpsd(var.getXmm(), var.getXmm(), int(asmjit::x86::kCmpEQ));
			jc->cc->andnpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		}
		return var;
	}

	double MpOperationIsNan::evaluateDRetD(double * args) { return std::isnan(args[0]) ? 1.0 : 0.0; }

	std::complex<double> MpOperationIsNan::evaluateCRetC(std::complex<double>* args)
	{
		return std::complex<double>(std::isnan(args[0].real()) ? 1.0 : 0.0, std::isnan(args[0].imag()) ? 1.0 : 0.0);
	}

	// MpOperationGetReal
	JitVar MpOperationGetReal::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar varRet(jc->cc->newXmmSd(), JitVar::FLAGS::FLAG_NONE);;

		if (node->takesComplex())
		{
			jc->cc->xorpd(varRet.getXmm(), varRet.getXmm());
			if (var.isXmm())
			{
				jc->cc->movsd(varRet.getXmm(), var.getXmm());
			}
			else
			{
				jc->cc->movsd(varRet.getXmm(), var.getMem());
			}
		}
		else
		{
			throw std::runtime_error("should not be reached");
		}
		return varRet;
	}

	double MpOperationGetReal::evaluateCRetD(std::complex<double>* args) { return args->real(); }

	// MpOperationGetImag
	JitVar MpOperationGetImag::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->onNode(node->getAt(0)));

		JitVar varRet(jc->cc->newXmmSd(), JitVar::FLAGS::FLAG_NONE);;
		if (node->takesComplex())
		{
			var = jc->registerVarComplex(var, !node->getAt(0)->hasNodeFlag(kAstReturnsComplex));
			jc->cc->xorpd(varRet.getXmm(), varRet.getXmm());
			jc->cc->shufpd(var.getXmm(), varRet.getXmm(), asmjit::x86::shufImm(0, 1));
		}
		else
		{
			throw std::runtime_error("should not be reached");
		}
		return var;
	}

	double MpOperationGetImag::evaluateCRetD(std::complex<double>* args) { return args->imag(); }

	// Square root
	JitVar mathpresso::MpOperationSqrt::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var = jc->onNode(node->getAt(0));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAGS::FLAG_NONE);
		if (var.isXmm())
			jc->cc->sqrtsd(result.getXmm(), var.getXmm());
		else
			jc->cc->sqrtsd(result.getXmm(), var.getMem());
		return result;
	}

	double MpOperationSqrt::evaluateDRetD(double * args) { return std::sqrt(args[0]); }

	// Square root, complex result
	std::complex<double> sqrtRC(double  * x) { return std::sqrt(std::complex<double>(x[0], 0)); }
	std::complex<double> sqrtCC(std::complex<double> *  x) { return std::sqrt(x[0]); }

	MpOperationSqrtC::MpOperationSqrtC() : MpOperationFunc(1, MpOperationFlags::OpFlagDReturnsC, reinterpret_cast<void*>(sqrtRC), reinterpret_cast<void*>(sqrtCC))
	{}

	// Negation
	JitVar MpOperationNeg::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var = jc->onNode(node->getAt(0));

		if (!hasFlag(kAstReturnsComplex))
		{
			var = jc->writableVar(var);
			jc->cc->pxor(var.getXmm(), jc->getConstantU64(uint64_t(0x8000000000000000)).getMem());
		}
		else
		{
			var = jc->writableVarComplex(var);
			jc->cc->pxor(var.getXmm(), jc->getConstantU64(uint64_t(0x8000000000000000), uint64_t(0x8000000000000000)).getMem());
		}
		return var;
	}

	double MpOperationNeg::evaluateDRetD(double * args) { return -args[0]; }

	std::complex<double> MpOperationNeg::evaluateCRetC(std::complex<double>* args) { return -args[0]; }

	uint32_t MpOperationNeg::optimizeSpecial(AstOptimizer * opt, AstNode * node)
	{
		// -(-(x)) = x
		if (node->getAt(0)->getNodeType() == kAstNodeUnaryOp
			&& static_cast<AstUnaryOp*>(node)->mpOp_ == static_cast<AstUnaryOp*>(node->getAt(0))->mpOp_)
		{
			AstNode* childOfChild = static_cast<AstUnaryOp*>(node->getAt(0))->unlinkChild();
			node->getParent()->replaceNode(node, childOfChild);
			opt->getAst()->deleteNode(node);

		}
		return kErrorOk;
	}

	// Not
	JitVar MpOperationNot::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var = jc->onNode(node->getAt(0));

		if (hasFlag(kAstReturnsComplex))
		{
			var = jc->writableVarComplex(var);
			jc->cc->cmppd(var.getXmm(), jc->getConstantD64(std::complex<double>(0.0, 0.0)).getMem(), int(asmjit::x86::kCmpEQ));
			jc->cc->andpd(var.getXmm(), jc->getConstantD64(std::complex<double>(1.0, 0.0)).getMem());
		}
		else
		{
			var = jc->writableVar(var);
			jc->cc->cmpsd(var.getXmm(), jc->getConstantD64AsPD(0.0).getMem(), int(asmjit::x86::kCmpEQ));
			jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		}
		return var;
	}

	double MpOperationNot::evaluateDRetD(double * args) { return args[0] == 0 ? 1.0 : 0.0; }

	std::complex<double> MpOperationNot::evaluateCRetC(std::complex<double>* args) { return std::complex<double>(args[0] == std::complex<double>(0, 0) ? 1.0 : 0.0, 0.0); }

	// Conjugate
	JitVar MpOperationConjug::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar tmp = jc->onNode(node->getAt(0));
		JitVar result = jc->registerVarComplex(tmp, !node->getAt(0)->returnsComplex());
		jc->cc->pxor(result.getXmm(), jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000)).getMem());
		return result;
	}

	std::complex<double> MpOperationConjug::evaluateCRetC(std::complex<double>* args) { return std::complex<double>(args->real(), -args->imag()); }

	uint32_t MpOperationConjug::optimizeSpecial(AstOptimizer * opt, AstNode * node)
	{
		// conj(conj(x)) = x
		if (node->getAt(0)->getNodeType() == kAstNodeUnaryOp
			&& static_cast<AstUnaryOp*>(node)->mpOp_ == static_cast<AstUnaryOp*>(node->getAt(0))->mpOp_)
		{
			AstNode* childOfChild = static_cast<AstUnaryOp*>(node->getAt(0))->unlinkChild();
			node->getParent()->replaceNode(node, childOfChild);
			opt->getAst()->deleteNode(node);
		}
		return kErrorOk;
	}

	// Reciprocal
	JitVar mathpresso::MpOperationRecip::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var = jc->onNode(node->getAt(0));
		JitVar result;
		if (node->takesComplex())
		{
			// as of http://www.chemistrylearning.com/reciprocal-of-a-complex-number/
			var = jc->writableVarComplex(var);
			result = JitVar(jc->cc->newXmmPd(), JitVar::FLAG_NONE);
			jc->cc->movapd(result.getXmm(), var.getXmm());
			jc->cc->mulpd(var.getXmm(), var.getXmm());
			jc->cc->haddpd(var.getXmm(), var.getXmm());
			jc->cc->pxor(result.getXmm(), jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000)).getMem());
			jc->cc->divpd(result.getXmm(), var.getXmm());
		}
		else
		{
			result = JitVar(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			jc->cc->movsd(result.getXmm(), jc->getConstantD64(1.0).getMem());
			if (var.isMem())
				jc->cc->divsd(result.getXmm(), var.getMem());
			else
				jc->cc->divsd(result.getXmm(), var.getXmm());
		}
		return result;
	}

	std::complex<double> MpOperationRecip::evaluateCRetC(std::complex<double>* args) { return 1.0 / args[0]; }

	double MpOperationRecip::evaluateDRetD(double * args) { return 1.0 / args[0]; }

	// MpOperationTrigonometrie
	// helpers:
	double _sin(double * arg) { return std::sin(arg[0]); }
	double _cos(double * arg) { return std::cos(arg[0]); }
	double _tan(double * arg) { return std::tan(arg[0]); }
	double _asin(double * arg) { return std::asin(arg[0]); }
	double _acos(double * arg) { return std::acos(arg[0]); }
	double _atan(double * arg) { return std::atan(arg[0]); }
	double _sinh(double * arg) { return std::sinh(arg[0]); }
	double _cosh(double * arg) { return std::cosh(arg[0]); }
	double _tanh(double * arg) { return std::tanh(arg[0]); }
	std::complex<double> _sinC(std::complex<double>* arg) { return std::sin(arg[0]); }
	std::complex<double> _cosC(std::complex<double>* arg) { return std::cos(arg[0]); }
	std::complex<double> _tanC(std::complex<double>* arg) { return std::tan(arg[0]); }
	std::complex<double> _asinC(std::complex<double>* arg) { return std::asin(arg[0]); }
	std::complex<double> _acosC(std::complex<double>* arg) { return std::acos(arg[0]); }
	std::complex<double> _atanC(std::complex<double>* arg) { return std::atan(arg[0]); }
	std::complex<double> _sinhC(std::complex<double>* arg) { return std::sinh(arg[0]); }
	std::complex<double> _coshC(std::complex<double>* arg) { return std::cosh(arg[0]); }
	std::complex<double> _tanhC(std::complex<double>* arg) { return std::tanh(arg[0]); }

	JitVar MpOperationTrigonometrie::compile(JitCompiler * jc, AstNode * node)
	{
		void* tmpC = fnC_;
		void* tmpD = fnD_;
		switch (type_)
		{
		case trigonometrieFunc::sin:
			fnD_ = reinterpret_cast<void*>(_sin);
			fnC_ = reinterpret_cast<void*>(_sinC);
			break;
		case trigonometrieFunc::cos:
			fnD_ = reinterpret_cast<void*>(_cos);
			fnC_ = reinterpret_cast<void*>(_cosC);
			break;
		case trigonometrieFunc::tan:
			fnD_ = reinterpret_cast<void*>(_tan);
			fnC_ = reinterpret_cast<void*>(_tanC);
			break;
		case trigonometrieFunc::asin:
			fnD_ = reinterpret_cast<void*>(_asin);
			fnC_ = reinterpret_cast<void*>(_asinC);
			break;
		case trigonometrieFunc::acos:
			fnD_ = reinterpret_cast<void*>(_acos);
			fnC_ = reinterpret_cast<void*>(_acosC);
			break;
		case trigonometrieFunc::atan:
			fnD_ = reinterpret_cast<void*>(_atan);
			fnC_ = reinterpret_cast<void*>(_atanC);
			break;
		case trigonometrieFunc::sinh:
			fnD_ = reinterpret_cast<void*>(_sinh);
			fnC_ = reinterpret_cast<void*>(_sinhC);
			break;
		case trigonometrieFunc::cosh:
			fnD_ = reinterpret_cast<void*>(_cosh);
			fnC_ = reinterpret_cast<void*>(_coshC);
			break;
		case trigonometrieFunc::tanh:
			fnD_ = reinterpret_cast<void*>(_tanh);
			fnC_ = reinterpret_cast<void*>(_tanhC);
			break;
		default:
			throw std::runtime_error("no function of this type available.");
		}
		JitVar ret = MpOperationFunc::compile(jc, node);
		fnD_ = tmpD;
		fnC_ = tmpC;
		return ret;
	}

	double MpOperationTrigonometrie::evaluateDRetD(double * args)
	{
		switch (type_)
		{
		case trigonometrieFunc::sin: return _sin(args);
		case trigonometrieFunc::cos: return _cos(args);
		case trigonometrieFunc::tan: return _tan(args);
		case trigonometrieFunc::asin: return _asin(args);
		case trigonometrieFunc::acos: return _acos(args);
		case trigonometrieFunc::atan: return _atan(args);
		case trigonometrieFunc::sinh: return _sinh(args);
		case trigonometrieFunc::cosh: return _cosh(args);
		case trigonometrieFunc::tanh: return _tanh(args);
		default:
			throw std::runtime_error("no function of this type available.");
		}
	}

	std::complex<double> MpOperationTrigonometrie::evaluateCRetC(std::complex<double>* args)
	{
		switch (type_)
		{
		case trigonometrieFunc::sin: return _sinC(args);
		case trigonometrieFunc::cos: return _cosC(args);
		case trigonometrieFunc::tan: return _tanC(args);
		case trigonometrieFunc::asin: return _asinC(args);
		case trigonometrieFunc::acos: return _acosC(args);
		case trigonometrieFunc::atan: return _atanC(args);
		case trigonometrieFunc::sinh: return _sinhC(args);
		case trigonometrieFunc::cosh: return _coshC(args);
		case trigonometrieFunc::tanh: return _tanhC(args);
		default:
			throw std::runtime_error("no function of this type available.");
		}
	}

	// sign bit
	JitVar MpOperationSignBit::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
		jc->cc->pshufd(result.getXmm(), jc->registerVar(var).getXmm(), asmjit::x86::shufImm(3, 2, 1, 1));
		jc->cc->psrad(result.getXmm(), 31);
		jc->cc->andpd(result.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return result;
	}

	double MpOperationSignBit::evaluateDRetD(double * args) { return std::signbit(args[0]) ? 1.0 : 0.0; }

	// Copy sign
	JitVar MpOperationCopySign::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar vl = jc->writableVar(jc->onNode(node->getAt(0)));
		JitVar vr = jc->writableVar(jc->onNode(node->getAt(1)));

		jc->cc->andpd(vl.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x7FFFFFFFFFFFFFFF)).getMem());
		jc->cc->andpd(vr.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->orpd(vl.getXmm(), vr.getXmm());

		return vl;
	}

	double MpOperationCopySign::evaluateDRetD(double * args) { return std::copysign(args[0], args[1]); }

	// Average
	JitVar MpOperationAvg::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar vl = jc->onNode(node->getAt(0));;
		JitVar vr = jc->onNode(node->getAt(1));

		if (node->takesComplex())
		{
			if (!node->getAt(0)->returnsComplex())
				vl = jc->registerVarAsComplex(vl);
			else
				vl = jc->writableVarComplex(vl);

			if (!node->getAt(1)->returnsComplex())
				vr = jc->registerVarAsComplex(vr);

			if (vr.isMem())
				jc->cc->addpd(vl.getXmm(), vr.getMem());
			else
				jc->cc->addpd(vl.getXmm(), vr.getXmm());
			jc->cc->mulpd(vl.getXmm(), jc->getConstantD64(std::complex<double>(0.5, 0.5)).getXmm());
		}
		else
		{
			vl = jc->writableVar(jc->onNode(node->getAt(0)));
			vr = jc->onNode(node->getAt(1));
			if (vr.isMem())
				jc->cc->addsd(vl.getXmm(), vr.getMem());
			else
				jc->cc->addsd(vl.getXmm(), vr.getXmm());
			jc->cc->mulsd(vl.getXmm(), jc->getConstantD64(0.5).getXmm());
		}
		return vl;
	}

	double MpOperationAvg::evaluateDRetD(double * args) { return (args[0] + args[1]) * 0.5; }

	std::complex<double> MpOperationAvg::evaluateCRetC(std::complex<double>* args) { return (args[0] + args[1]) * 0.5; }

	// Absolute
	double absc(std::complex<double>* arg) { return std::abs(arg[0]); }
	JitVar MpOperationAbs::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar result;
		if (node->takesComplex())
		{
			fnC_ = reinterpret_cast<void*>(absc);
			result = MpOperationFunc::compile(jc, node);
		}
		else
		{
			var = jc->writableVar(var);
			result = JitVar(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			jc->cc->xorpd(result.getXmm(), result.getXmm());
			jc->cc->subsd(result.getXmm(), var.getXmm());
			jc->cc->maxsd(result.getXmm(), var.getXmm());
		}
		return result;
	}

	double MpOperationAbs::evaluateDRetD(double * args) { return std::abs(args[0]); }

	double MpOperationAbs::evaluateCRetD(std::complex<double>* args) { return std::abs(args[0]); }

	// round
	JitVar MpOperationRound::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

		if (jc->enableSSE4_1)
		{
			JitVar tmp(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			jc->cc->roundsd(tmp.getXmm(), var.getXmm(), asmjit::x86::kRoundDown | asmjit::x86::kRoundInexact);
			if (result.getXmm().getId() != var.getXmm().getId())
				jc->cc->movsd(result.getXmm(), var.getXmm());
			jc->cc->subsd(result.getXmm(), tmp.getXmm());
			jc->cc->cmpsd(result.getXmm(), jc->getConstantD64(0.5).getMem(), asmjit::x86::kCmpNLT);
			jc->cc->andpd(result.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
			jc->cc->addpd(result.getXmm(), tmp.getXmm());
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;
			const double magic1 = 6755399441055745.0;

			JitVar t1(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t2(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t3(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

			jc->cc->movsd(t2.getXmm(), var.getXmm());
			jc->cc->movsd(t3.getXmm(), var.getXmm());
			if (result.getXmm().getId() != var.getXmm().getId())
				jc->cc->movsd(result.getXmm(), var.getXmm());
			jc->cc->addsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->addsd(t3.getXmm(), jc->getConstantD64(magic1).getMem());
			jc->cc->movsd(t1.getXmm(), var.getXmm());
			jc->cc->subsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->subsd(t3.getXmm(), jc->getConstantD64(magic1).getMem());
			jc->cc->cmpsd(t1.getXmm(), jc->getConstantD64(maxn).getMem(), asmjit::x86::kCmpNLT);
			jc->cc->maxsd(t2.getXmm(), t3.getXmm());
			jc->cc->andpd(result.getXmm(), t1.getXmm());
			jc->cc->andnpd(t1.getXmm(), t2.getXmm());
			jc->cc->orpd(result.getXmm(), t1.getXmm());
		}

		return result;
	}

	double MpOperationRound::evaluateDRetD(double * args)
	{
		double y = ::floor(args[0]);
		return y + (args[0] - y >= 0.5 ? double(1.0) : double(0.0));
	}

	// roundeven
	JitVar MpOperationRoundEven::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(result.getXmm(), var.getXmm(), asmjit::x86::kRoundNearest | asmjit::x86::kRoundInexact);
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t2(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

			jc->cc->movsd(t1.getXmm(), var.getXmm());
			jc->cc->movsd(t2.getXmm(), var.getXmm());
			jc->cc->addsd(t1.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->cmpsd(t2.getXmm(), jc->getConstantD64(maxn).getMem(), asmjit::x86::kCmpNLT);
			jc->cc->subsd(t1.getXmm(), jc->getConstantD64(magic0).getMem());
			if (result.getXmm().getId() != var.getXmm().getId())
				jc->cc->movsd(result.getXmm(), var.getXmm());
			jc->cc->andpd(result.getXmm(), t2.getXmm());
			jc->cc->andnpd(t2.getXmm(), t1.getXmm());
			jc->cc->orpd(result.getXmm(), t2.getXmm());
		}
		return result;
	}

	double MpOperationRoundEven::evaluateDRetD(double * args) { return std::rint(args[0]); }


	// trunc
	JitVar MpOperationTrunc::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(result.getXmm(), var.getXmm(), asmjit::x86::kRoundTrunc | asmjit::x86::kRoundInexact);
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t2(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t3(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

			jc->cc->movsd(t2.getXmm(), jc->getConstantU64(ASMJIT_UINT64_C(0x7FFFFFFFFFFFFFFF)).getMem());
			jc->cc->andpd(t2.getXmm(), var.getXmm());
			if (result.getXmm().getId() != var.getXmm().getId())
				jc->cc->movsd(result.getXmm(), var.getXmm());
			jc->cc->movsd(t1.getXmm(), t2.getXmm());
			jc->cc->addsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->movsd(t3.getXmm(), t1.getXmm());
			jc->cc->subsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->cmpsd(t1.getXmm(), jc->getConstantD64(maxn).getMem(), asmjit::x86::kCmpNLT);
			jc->cc->cmpsd(t3.getXmm(), t2.getXmm(), asmjit::x86::kCmpLT);
			jc->cc->orpd(t1.getXmm(), jc->getConstantU64AsPD(ASMJIT_UINT64_C(0x8000000000000000)).getMem());
			jc->cc->andpd(t3.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
			jc->cc->andpd(result.getXmm(), t1.getXmm());
			jc->cc->subpd(t2.getXmm(), t3.getXmm());
			jc->cc->andnpd(t1.getXmm(), t2.getXmm());
			jc->cc->orpd(result.getXmm(), t1.getXmm());
		}
		return result;
	}

	double MpOperationTrunc::evaluateDRetD(double * args) { return std::trunc(args[0]); }

	// frac
	JitVar MpOperationFrac::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar tmp(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(tmp.getXmm(), var.getXmm(), int(asmjit::x86::kRoundDown | asmjit::x86::kRoundInexact));
			jc->cc->subsd(var.getXmm(), tmp.getXmm());
			return var;
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t2(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t3(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

			jc->cc->movsd(t2.getXmm(), var.getXmm());
			jc->cc->movsd(t3.getXmm(), var.getXmm());
			if (tmp.getXmm().getId() != var.getXmm().getId())
				jc->cc->movsd(tmp.getXmm(), var.getXmm());
			jc->cc->addsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->movsd(t1.getXmm(), var.getXmm());
			jc->cc->subsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->cmpsd(t1.getXmm(), jc->getConstantD64(maxn).getMem(), asmjit::x86::kCmpNLT);
			jc->cc->cmpsd(t3.getXmm(), t2.getXmm(), asmjit::x86::kCmpLT);
			jc->cc->andpd(t3.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
			jc->cc->andpd(tmp.getXmm(), t1.getXmm());
			jc->cc->subpd(t2.getXmm(), t3.getXmm());
			jc->cc->andnpd(t1.getXmm(), t2.getXmm());
			jc->cc->orpd(tmp.getXmm(), t1.getXmm());

			jc->cc->subsd(var.getXmm(), tmp.getXmm());
		}
		return var;
	}

	double MpOperationFrac::evaluateDRetD(double * args) { return args[0] - std::floor(args[0]); }

	// floor
	JitVar MpOperationFloor::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(result.getXmm(), var.getXmm(), asmjit::x86::kRoundDown | asmjit::x86::kRoundInexact);
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t2(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t3(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

			jc->cc->movsd(t2.getXmm(), var.getXmm());
			jc->cc->movsd(t3.getXmm(), var.getXmm());
			if (result.getXmm().getId() != var.getXmm().getId())
				jc->cc->movsd(result.getXmm(), var.getXmm());
			jc->cc->addsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->movsd(t1.getXmm(), var.getXmm());
			jc->cc->subsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->cmpsd(t1.getXmm(), jc->getConstantD64(maxn).getMem(), asmjit::x86::kCmpNLT);
			jc->cc->cmpsd(t3.getXmm(), t2.getXmm(), asmjit::x86::kCmpLT);
			jc->cc->andpd(t3.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
			jc->cc->andpd(result.getXmm(), t1.getXmm());
			jc->cc->subpd(t2.getXmm(), t3.getXmm());
			jc->cc->andnpd(t1.getXmm(), t2.getXmm());
			jc->cc->orpd(result.getXmm(), t1.getXmm());
		}

		return result;
	}

	double MpOperationFloor::evaluateDRetD(double * args) { return std::floor(args[0]); }

	// ceil
	JitVar MpOperationcCeil::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(result.getXmm(), var.getXmm(), asmjit::x86::kRoundUp | asmjit::x86::kRoundInexact);
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t2(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t3(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

			jc->cc->movsd(t2.getXmm(), var.getXmm());
			jc->cc->movsd(t3.getXmm(), var.getXmm());
			if (result.getXmm().getId() != var.getXmm().getId())
				jc->cc->movsd(result.getXmm(), var.getXmm());
			jc->cc->addsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->movsd(t1.getXmm(), var.getXmm());
			jc->cc->subsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->cmpsd(t1.getXmm(), jc->getConstantD64(maxn).getMem(), asmjit::x86::kCmpNLT);
			jc->cc->cmpsd(t3.getXmm(), t2.getXmm(), asmjit::x86::kCmpNLE);
			jc->cc->andpd(t3.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
			jc->cc->andpd(result.getXmm(), t1.getXmm());
			jc->cc->addpd(t2.getXmm(), t3.getXmm());
			jc->cc->andnpd(t1.getXmm(), t2.getXmm());
			jc->cc->orpd(result.getXmm(), t1.getXmm());
		}

		return result;
	}

	double MpOperationcCeil::evaluateDRetD(double * args) { return std::ceil(args[0]); }

	// Log
	double logRR(double  * x) { return std::log(x[0]); }
	std::complex<double> logCC(std::complex<double> *  x) { return std::log(x[0]); }

	MpOperationLog::MpOperationLog() : MpOperationFunc(1, MpOperationFlags::OpFlagNone, reinterpret_cast<void*>(logRR), reinterpret_cast<void*>(logCC))
	{}

	// Log2
	double log2RR(double  * x) { return std::log2(x[0]); }
	std::complex<double> log2CC(std::complex<double> *  x) { return std::log(x[0]) / log(2); }

	MpOperationLog2::MpOperationLog2() : MpOperationFunc(1, MpOperationFlags::OpFlagNone, reinterpret_cast<void*>(log2RR), reinterpret_cast<void*>(log2CC))
	{}

	// Log10
	double log10RR(double  * x) { return std::log10(x[0]); }
	std::complex<double> log10CC(std::complex<double> *  x) { return std::log10(x[0]); }

	MpOperationLog10::MpOperationLog10() : MpOperationFunc(1, MpOperationFlags::OpFlagNone, reinterpret_cast<void*>(log10RR), reinterpret_cast<void*>(log10CC))
	{}

	// exp
	double expRR(double  * x) { return std::exp(x[0]); }
	std::complex<double> expCC(std::complex<double> *  x) { return std::exp(x[0]); }

	MpOperationExp::MpOperationExp() : MpOperationFunc(1, MpOperationFlags::OpFlagNone, reinterpret_cast<void*>(expRR), reinterpret_cast<void*>(expCC))
	{}

	// pow
	double powRR(double * x) { return std::pow(x[0], x[1]); }
	std::complex<double> powCC(std::complex<double> *  x) { return std::pow(x[0], x[1]); }

	MpOperationPow::MpOperationPow() : MpOperationFunc(2, MpOperationFlags::OpFlagNopIfROne, reinterpret_cast<void*>(powRR), reinterpret_cast<void*>(powCC))
	{}

	// Atan2
	double atan2RR(double * x) { return std::atan2(x[0], x[1]); }

	MpOperationAtan2::MpOperationAtan2() : MpOperationFunc(2, MpOperationFlags::OpFlagNone, reinterpret_cast<void*>(atan2RR), nullptr)
	{}

	// hypot
	double hypotRR(double * x) { return std::hypot(x[0], x[1]); }

	MpOperationHypot::MpOperationHypot() : MpOperationFunc(2, MpOperationFlags::OpFlagNone, reinterpret_cast<void*>(hypotRR), nullptr)
	{}

	// mpOperationBinary
	JitVar MpOperationBinary::compile(JitCompiler* jc, AstNode * node)
	{
		JitVar vl, vr;
		AstNode* left = node->getAt(0);
		AstNode* right = node->getAt(1);

		if (node->takesComplex())
		{
			// check whether the vars are the same, to reduce memory-operations
			if (left->isVar() && right->isVar() &&
				static_cast<AstVar*>(left)->getSymbol() == static_cast<AstVar*>(right)->getSymbol())
			{
				vl = vr = jc->writableVarComplex(jc->onNode(left));
			}
			else
			{
				// check that every node has onNode called on it and that they are registered as complex.
				// also make sure, vl is in a register.
				if (!left->returnsComplex())
					vl = jc->registerVarAsComplex(jc->onNode(left));
				else
					vl = jc->writableVarComplex(jc->onNode(left));

				if (!right->returnsComplex())
					vr = jc->registerVarAsComplex(jc->onNode(right));
				else
					vr = jc->onNode(right);
			}
			return generateAsmComplex(jc, vl, vr);
		}
		else
		{
			// check whether the vars are the same, to reduce memory-operations
			if (left->isVar() && right->isVar() &&
				static_cast<AstVar*>(left)->getSymbol() == static_cast<AstVar*>(right)->getSymbol())
			{
				vl = vr = jc->writableVar(jc->onNode(left));
			}
			else
			{
				// check that every node has onNode called on it and vl is in a register
				vl = jc->writableVar(jc->onNode(left));
				vr = jc->onNode(right);
			}
			return generatAsmReal(jc, vl, vr);
		}
	}

	uint32_t MpOperationBinary::optimize(AstOptimizer *opt, AstNode *node)
	{
		MATHPRESSO_PROPAGATE(opt->onNode(node->getAt(0)));
		AstNode* left = node->getAt(0);

		MATHPRESSO_PROPAGATE(opt->onNode(node->getAt(1)));
		AstNode* right = node->getAt(1);

		bool lIsImm = left->isImm();
		bool rIsImm = right->isImm();
		bool needs_complex = left->returnsComplex() || right->returnsComplex();

		// set the flags according to the operands and the capability's of the MpOperationbinary.
		if ((needs_complex || hasFlag(MpOperationFlags::OpHasNoReal)) && !hasFlag(MpOperationFlags::OpHasNoComplex))
		{
			node->addNodeFlags(AstNodeFlags::kAstReturnsComplex | AstNodeFlags::kAstTakesComplex);
		}

		if (lIsImm && rIsImm)
		{
			// optimize a calculation with two immediates.
			AstImm* lNode = static_cast<AstImm*>(left);
			AstImm* rNode = static_cast<AstImm*>(right);
			if (needs_complex && !hasFlag(MpOperationFlags::OpHasNoComplex))
			{
				lNode->setValue(calculateComplex(lNode->getValueCplx(), rNode->getValueCplx()));
			}
			else if (!needs_complex && !hasFlag(MpOperationFlags::OpHasNoReal))
			{
				lNode->setValue(calculateReal(lNode->getValue(), rNode->getValue()));
			}
			else if (!needs_complex && !hasFlag(MpOperationFlags::OpHasNoComplex))
			{
				lNode->setValue(calculateComplex(lNode->getValueCplx(), rNode->getValueCplx()));
			}
			else
			{
				throw std::runtime_error("Wrong implementation of MpOperationBinary!");
			}

			// setValue sets the correct flags automatically.
			node->_children[0]->_parent = nullptr;
			node->_children[0] = nullptr;
			node->getParent()->replaceNode(node, lNode);
			opt->getAst()->deleteNode(node);
		}
		else if (lIsImm)
		{
			AstImm* lNode = static_cast<AstImm*>(left);
			// if the node is real, the imaginary part is set to zero by default.
			if ((hasFlag(MpOperationFlags::OpFlagNopIfLZero) && lNode->getValueCplx() == std::complex<double>(0.0, 0.0)) ||
				(hasFlag(MpOperationFlags::OpFlagNopIfLOne) && lNode->getValueCplx() == std::complex<double>(1.0, 0.0)))
			{
				node->_children[1]->_parent = nullptr;
				node->_children[1] = nullptr;
				node->getParent()->replaceNode(node, right);
				opt->getAst()->deleteNode(node);
			}
		}
		else if (rIsImm)
		{
			AstImm* rNode = static_cast<AstImm*>(right);

			if ((hasFlag(MpOperationFlags::OpFlagNopIfRZero) && rNode->getValueCplx() == std::complex<double>(0.0, 0.0)) ||
				(hasFlag(MpOperationFlags::OpFlagNopIfROne) && rNode->getValueCplx() == std::complex<double>(1.0, 0.0)))
			{
				node->_children[0]->_parent = nullptr;
				node->_children[0] = nullptr;
				node->getParent()->replaceNode(node, left);
				opt->getAst()->deleteNode(node);
			}

		}
		return ErrorCode::kErrorOk;
	}

	JitVar MpOperationBinary::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) { throw std::runtime_error("No Override available!"); }

	JitVar MpOperationBinary::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) { throw std::runtime_error("No Override available!"); }

	// Addition
	JitVar MpOperationAdd::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->addsd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->addsd(vl.getXmm(), vr.getXmm());
		}
		return vl;

	}

	JitVar MpOperationAdd::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.isMem())
		{
			jc->cc->addpd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->addpd(vl.getXmm(), vr.getXmm());
		}
		return vl;
	}

	double MpOperationAdd::calculateReal(double vl, double vr) { return vl + vr; }
	std::complex<double> MpOperationAdd::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return vl + vr; }

	// Subtraction
	JitVar MpOperationSub::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->subsd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->subsd(vl.getXmm(), vr.getXmm());
		}
		return vl;

	}

	JitVar MpOperationSub::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.isMem())
		{
			jc->cc->subpd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->subpd(vl.getXmm(), vr.getXmm());
		}
		return vl;
	}

	double MpOperationSub::calculateReal(double vl, double vr) { return vl - vr; }
	std::complex<double> MpOperationSub::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return vl - vr; }

	// Multiplication
	JitVar MpOperationMul::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->mulsd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->mulsd(vl.getXmm(), vr.getXmm());
		}
		return vl;

	}

	JitVar MpOperationMul::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.isMem())
		{
			vr = jc->writableVarComplex(vr);
		}
		if (vl == vr)
		{
			vr = jc->copyVarComplex(vl, JitVar::FLAG_NONE);
		}
		JitVar ret(jc->cc->newXmmPd(), JitVar::FLAG_NONE);
		JitVar negateImag = jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000));

		// algorithm with modifications taken from:
		// https://www.codeproject.com/Articles/874396/Crunching-Numbers-with-AVX-and-AVX
		jc->cc->movapd(ret.getXmm(), vl.getXmm());
		jc->cc->mulpd(ret.getXmm(), vr.getXmm());
		jc->cc->shufpd(vr.getXmm(), vr.getXmm(), asmjit::x86::shufImm(0, 1));
		jc->cc->pxor(vr.getXmm(), negateImag.getMem());
		jc->cc->mulpd(vl.getXmm(), vr.getXmm());
		jc->cc->hsubpd(ret.getXmm(), vl.getXmm());
		return ret;
	}

	double MpOperationMul::calculateReal(double vl, double vr) { return vl * vr; }
	std::complex<double> MpOperationMul::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return vl * vr; }

	// Division
	JitVar MpOperationDiv::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->divsd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->divsd(vl.getXmm(), vr.getXmm());
		}
		return vl;

	}

	JitVar MpOperationDiv::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.isMem())
		{
			vr = jc->writableVarComplex(vr);
		}
		if (vl == vr)
		{
			vr = jc->copyVarComplex(vl, JitVar::FLAG_NONE);
		}
		JitVar ret(jc->cc->newXmmPd(), JitVar::FLAG_NONE);
		JitVar negateImag = jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000));

		jc->cc->pxor(vr.getXmm(), negateImag.getMem());

		jc->cc->movapd(ret.getXmm(), vl.getXmm());
		jc->cc->mulpd(ret.getXmm(), vr.getXmm());
		jc->cc->shufpd(vr.getXmm(), vr.getXmm(), asmjit::x86::shufImm(0, 1));
		jc->cc->pxor(vr.getXmm(), negateImag.getMem());
		jc->cc->mulpd(vl.getXmm(), vr.getXmm());
		jc->cc->hsubpd(ret.getXmm(), vl.getXmm());

		jc->cc->mulpd(vr.getXmm(), vr.getXmm());
		jc->cc->haddpd(vr.getXmm(), vr.getXmm());
		jc->cc->divpd(ret.getXmm(), vr.getXmm());
		return ret;
	}

	double MpOperationDiv::calculateReal(double vl, double vr) { return vl / vr; }
	std::complex<double> MpOperationDiv::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return vl / vr; }

	// Minimum
	JitVar MpOperationMin::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->minsd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->minsd(vl.getXmm(), vr.getXmm());
		}
		return vl;
	}

	double MpOperationMin::calculateReal(double vl, double vr) { return std::min(vl, vr); }

	// Maximum
	JitVar MpOperationMax::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->maxsd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->maxsd(vl.getXmm(), vr.getXmm());
		}
		return vl;
	}

	double MpOperationMax::calculateReal(double vl, double vr) { return std::max(vl, vr); }

	// Equality
	JitVar MpOperationEq::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getMem(), asmjit::x86::kCmpEQ);
		}
		else
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getXmm(), asmjit::x86::kCmpEQ);
		}
		jc->cc->andpd(vl.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return vl;

	}

	JitVar MpOperationEq::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.isMem())
		{
			jc->cc->cmppd(vl.getXmm(), vr.getMem(), asmjit::x86::kCmpEQ);
		}
		else
		{
			jc->cc->cmppd(vl.getXmm(), vr.getXmm(), asmjit::x86::kCmpEQ);
		}
		jc->cc->haddpd(vl.getXmm(), vl.getXmm());
		jc->cc->andpd(vl.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return vl;
	}

	double MpOperationEq::calculateReal(double vl, double vr) { return vl == vr ? 1.0 : 0.0; }
	std::complex<double> MpOperationEq::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return std::complex<double>(vl == vr ? 1.0 : 0.0, 0.0); }

	// Inequality
	JitVar MpOperationNe::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getMem(), asmjit::x86::kCmpNEQ);
		}
		else
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getXmm(), asmjit::x86::kCmpNEQ);
		}
		jc->cc->andpd(vl.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return vl;

	}

	JitVar MpOperationNe::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.isMem())
		{
			jc->cc->cmppd(vl.getXmm(), vr.getMem(), asmjit::x86::kCmpNEQ);
		}
		else
		{
			jc->cc->cmppd(vl.getXmm(), vr.getXmm(), asmjit::x86::kCmpNEQ);
		}
		jc->cc->haddpd(vl.getXmm(), vl.getXmm());
		jc->cc->andpd(vl.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return vl;
	}

	double MpOperationNe::calculateReal(double vl, double vr) { return vl != vr ? 1.0 : 0.0; }
	std::complex<double> MpOperationNe::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return std::complex<double>(vl != vr ? 1.0 : 0.0, 0.0); }

	// Lesser than
	JitVar MpOperationLt::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getMem(), asmjit::x86::kCmpLT);
		}
		else
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getXmm(), asmjit::x86::kCmpLT);
		}
		jc->cc->andpd(vl.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return vl;

	}

	double MpOperationLt::calculateReal(double vl, double vr) { return vl < vr ? 1.0 : 0.0; }

	// Lesser equal
	JitVar MpOperationLe::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getMem(), asmjit::x86::kCmpLE);
		}
		else
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getXmm(), asmjit::x86::kCmpLE);
		}
		jc->cc->andpd(vl.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return vl;

	}

	double MpOperationLe::calculateReal(double vl, double vr) { return vl <= vr ? 1.0 : 0.0; }

	// Greater than
	JitVar MpOperationGt::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getMem(), asmjit::x86::kCmpNLE);
		}
		else
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getXmm(), asmjit::x86::kCmpNLE);
		}
		jc->cc->andpd(vl.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return vl;

	}

	double MpOperationGt::calculateReal(double vl, double vr) { return vl > vr ? 1.0 : 0.0; }

	// Greater equal
	JitVar MpOperationGe::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		if (vr.getOperand().isMem())
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getMem(), asmjit::x86::kCmpNLT);
		}
		else
		{
			jc->cc->cmpsd(vl.getXmm(), vr.getXmm(), asmjit::x86::kCmpNLT);
		}
		jc->cc->andpd(vl.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return vl;

	}

	double MpOperationGe::calculateReal(double vl, double vr) { return vl >= vr ? 1.0 : 0.0; }

	// Modulo
	JitVar MpOperationModulo::generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
		JitVar tmp(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

		vl = jc->writableVar(vl);
		if (vl == vr)
		{
			vr = jc->copyVar(vl, JitVar::FLAG_NONE);
		}
		else
		{
			vr = jc->registerVar(vr);
		}

		jc->cc->movsd(result.getXmm(), vl.getXmm());
		jc->cc->divsd(vl.getXmm(), vr.getXmm());

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(tmp.getXmm(), vl.getXmm(), asmjit::x86::kRoundTrunc | asmjit::x86::kRoundInexact);
		}
		else
		{
			JitVar var = vl;

			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t2(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
			JitVar t3(jc->cc->newXmmSd(), JitVar::FLAG_NONE);

			jc->cc->movsd(t2.getXmm(), jc->getConstantU64(ASMJIT_UINT64_C(0x7FFFFFFFFFFFFFFF)).getMem());
			jc->cc->andpd(t2.getXmm(), var.getXmm());
			if (result.getXmm().getId() != var.getXmm().getId())
				jc->cc->movsd(tmp.getXmm(), var.getXmm());
			jc->cc->movsd(t1.getXmm(), t2.getXmm());
			jc->cc->addsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->movsd(t3.getXmm(), t1.getXmm());
			jc->cc->subsd(t2.getXmm(), jc->getConstantD64(magic0).getMem());
			jc->cc->cmpsd(t1.getXmm(), jc->getConstantD64(maxn).getMem(), asmjit::x86::kCmpNLT);
			jc->cc->cmpsd(t3.getXmm(), t2.getXmm(), asmjit::x86::kCmpLT);
			jc->cc->orpd(t1.getXmm(), jc->getConstantU64AsPD(ASMJIT_UINT64_C(0x8000000000000000)).getMem());
			jc->cc->andpd(t3.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
			jc->cc->andpd(tmp.getXmm(), t1.getXmm());
			jc->cc->subpd(t2.getXmm(), t3.getXmm());
			jc->cc->andnpd(t1.getXmm(), t2.getXmm());
			jc->cc->orpd(tmp.getXmm(), t1.getXmm());


		}

		jc->cc->mulsd(tmp.getXmm(), vr.getXmm());
		jc->cc->subsd(result.getXmm(), tmp.getXmm());

		return result;
	}

	double MpOperationModulo::calculateReal(double vl, double vr) { return fmod(vl, vr); }

	// Ternary operation
	JitVar MpOperationTernary::compile(JitCompiler* jc, AstNode * node)
	{
		asmjit::Label lblElse = jc->cc->newLabel();
		asmjit::Label lblEnd = jc->cc->newLabel();
		JitVar erg;
		AstNode* left = static_cast<AstTernaryOp*>(node)->getLeft();
		AstNode* right = static_cast<AstTernaryOp*>(node)->getRight();
		AstNode* condition = static_cast<AstTernaryOp*>(node)->getCondition();

		JitVar ret = jc->onNode(condition);

		if (condition->returnsComplex())
			jc->cc->haddpd(ret.getXmm(), ret.getXmm());

		jc->cc->ucomisd(ret.getXmm(), jc->getConstantD64(0).getMem());
		jc->cc->je(lblElse);

		asmjit::X86Xmm regErg = jc->cc->newXmmPd();
		JitVar ergLeft = jc->onNode(left);

		bool lIsVarOrImm = left->getNodeType() == kAstNodeVar || left->getNodeType() == kAstNodeImm;
		bool rIsVarOrImm = right->getNodeType() == kAstNodeVar || right->getNodeType() == kAstNodeImm;

		if (lIsVarOrImm)
		{
			if (left->returnsComplex())
			{
				jc->cc->movupd(regErg, ergLeft.getXmm());
			}
			else
			{
				jc->cc->xorpd(regErg, regErg);
				jc->cc->movsd(regErg, ergLeft.getXmm());
			}
		}

		jc->cc->jmp(lblEnd);
		jc->cc->bind(lblElse);

		JitVar ergRight = jc->onNode(right);

		if (rIsVarOrImm)
		{
			if (right->returnsComplex())
			{
				jc->cc->movupd(regErg, ergRight.getXmm());
			}
			else
			{
				jc->cc->xorpd(regErg, regErg);
				jc->cc->movsd(regErg, ergRight.getXmm());
			}
		}

		jc->cc->bind(lblEnd);

		return jc->copyVarComplex(JitVar(regErg, JitVar::FLAG_NONE), JitVar::FLAG_NONE);
	}

	uint32_t MpOperationTernary::optimize(AstOptimizer *opt, AstNode *node)
	{
		AstBinaryOp* lastColon = static_cast<AstBinaryOp*>(node);
		// go to the last Colon after question-marks.
		while (lastColon->getOp() == kOpQMark)
		{
			lastColon = static_cast<AstBinaryOp*>(lastColon->getRight());
		}
		while (lastColon->getRight()->getOp() == kOpColon)
		{
			lastColon = static_cast<AstBinaryOp*>(lastColon->getRight());
		}

		if (lastColon->getOp() != kOpColon)
		{
			return opt->_errorReporter->onError(kErrorInvalidSyntax, node->getPosition(),
				"Invalid ternary operation. Expected a ':', found '%s' instead.", OpInfo::get(lastColon->getOp()).name);
		}

		AstNode* branchCondition = static_cast<AstBinaryOp*>(node)->getLeft();
		AstNode* branchLeft = lastColon->getLeft();
		AstNode* branchRight = lastColon->getRight();

		// remove branchCondition from the AST
		static_cast<AstBinaryOp*>(node)->setLeft(nullptr);
		branchCondition->_parent = nullptr;

		// remove the right path from the AST.
		lastColon->setRight(nullptr);
		branchRight->_parent = nullptr;


		// Distinguish between a complex and a non-complex case:
		// i.e.: cond1 ? cond2 ? a : b : c
		if (static_cast<AstBinaryOp*>(node)->getRight() != lastColon)
		{
			// remove left branch from the AST.
			branchLeft = static_cast<AstBinaryOp*>(node)->getRight();
			static_cast<AstBinaryOp*>(node)->setRight(nullptr);
			branchLeft->_parent = nullptr;

			// correct the right path.
			AstBinaryOp* preLastColon = static_cast<AstBinaryOp*>(lastColon->getParent());
			preLastColon->replaceAt(1, lastColon->getLeft());

		}
		// i.e.: cond1 ? a : b
		else
		{
			// remove left branch from the AST.
			lastColon->setLeft(nullptr);
			branchLeft->_parent = nullptr;
		}

		// create the new Ternary Node.
		AstTernaryOp* ternaryNode = node->getAst()->newNode<AstTernaryOp>(kOpQMark);
		ternaryNode->setCondition(branchCondition);
		ternaryNode->setLeft(branchLeft);
		ternaryNode->setRight(branchRight);
		ternaryNode->mpOp_ = opt->_symbols->at("?$2").get();

		AstBinaryOp* oldNode = static_cast<AstBinaryOp*>(node);

		// add the new node to the AST.
		node->getParent()->replaceNode(node, ternaryNode);

		// clean up:
		lastColon->setLeft(nullptr);
		opt->_ast->deleteNode(lastColon);
		opt->_ast->deleteNode(node);

		MATHPRESSO_PROPAGATE(opt->onNode(ternaryNode->getCondition()));
		AstNode* branchCond = ternaryNode->getCondition();
		if (branchCond->isImm())
		{
			// optimize an immediate condition
			bool conditionIsTrue = static_cast<AstImm*>(branchCond)->getValueCplx() != std::complex<double>({ 0, 0 });

			AstNode* nodeOptimized;

			if (conditionIsTrue)
			{
				nodeOptimized = ternaryNode->getLeft();
				ternaryNode->setLeft(nullptr);
			}
			else
			{
				nodeOptimized = ternaryNode->getRight();
				ternaryNode->setRight(nullptr);
			}

			nodeOptimized->_parent = nullptr;
			ternaryNode->getParent()->replaceNode(ternaryNode, nodeOptimized);
			ternaryNode->setCondition(nullptr);
			ternaryNode->setLeft(nullptr);
			ternaryNode->setRight(nullptr);
			opt->_ast->deleteNode(ternaryNode);
			MATHPRESSO_PROPAGATE(opt->onNode(nodeOptimized));

		}
		else
		{
			ternaryNode->removeNodeFlags(kAstTakesComplex | kAstReturnsComplex);
			MATHPRESSO_PROPAGATE(opt->onNode(ternaryNode->getLeft()));
			MATHPRESSO_PROPAGATE(opt->onNode(ternaryNode->getRight()));
			bool needs_complex = ternaryNode->getLeft()->returnsComplex() | ternaryNode->getRight()->returnsComplex();
			if (needs_complex)
			{
				ternaryNode->addNodeFlags(kAstReturnsComplex | kAstTakesComplex);
			}
		}
		return kErrorOk;

	}

	// Assignment
	JitVar mathpresso::MpOperationAssignment::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar result;
		AstVarDecl * varDecl = static_cast<AstVarDecl*>(node);

		if (varDecl->hasChild())
			result = jc->onNode(varDecl->getChild());

		AstSymbol* sym = varDecl->getSymbol();
		uint32_t slotId = sym->getVarSlotId();

		result.setRO();
		jc->varSlots[slotId] = result;

		return result;
	}

	uint32_t mathpresso::MpOperationAssignment::optimize(AstOptimizer * opt, AstNode * node)
	{
		AstVarDecl * varDecl;
		if (node->getNodeType() == kAstNodeVarDecl)
		{
			varDecl = static_cast<AstVarDecl*>(node);
		}
		else
		{
			return kErrorInvalidState;
		}

		AstSymbol* sym = varDecl->getSymbol();

		if (varDecl->hasChild())
		{
			MATHPRESSO_PROPAGATE(opt->onNode(varDecl->getChild()));
			AstNode* child = varDecl->getChild();

			if (child->returnsComplex())
			{
				varDecl->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);
				sym->setSymbolFlag(kAstSymbolIsComplex);
			}

			if (child->isImm())
			{
				if (child->returnsComplex())
				{
					sym->setValue(static_cast<AstImm*>(child)->getValueCplx());
				}
				else
				{
					sym->setValue(static_cast<AstImm*>(child)->getValue());
				}
				sym->setAssigned();
			}
		}

		return kErrorOk;
	}

} // end namespace mathpresso
