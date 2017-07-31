// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

#include <mathpresso/mpoperation_p.h>
#include <mathpresso/mpast_p.h>
#include <mathpresso/mpcompiler_p.h>
#include <mathpresso/mpoptimizer_p.h>
#include <asmjit/x86/x86operand.h>
#include <asmjit/x86/x86inst.h>

#include <complex>

namespace mathpresso {

#define VPTR(function) reinterpret_cast<void*>(function)

	// helpers:
	double _sin(double arg) { return std::sin(arg); }
	double _cos(double arg) { return std::cos(arg); }
	double _tan(double arg) { return std::tan(arg); }
	double _asin(double arg) { return std::asin(arg); }
	double _acos(double arg) { return std::acos(arg); }
	double _atan(double arg) { return std::atan(arg); }
	double _sinh(double arg) { return std::sinh(arg); }
	double _cosh(double arg) { return std::cosh(arg); }
	double _tanh(double arg) { return std::tanh(arg); }
	std::complex<double> _sinC(std::complex<double>* arg) { return std::sin(arg[0]); }
	std::complex<double> _cosC(std::complex<double>* arg) { return std::cos(arg[0]); }
	std::complex<double> _tanC(std::complex<double>* arg) { return std::tan(arg[0]); }
	std::complex<double> _asinC(std::complex<double>* arg) { return std::asin(arg[0]); }
	std::complex<double> _acosC(std::complex<double>* arg) { return std::acos(arg[0]); }
	std::complex<double> _atanC(std::complex<double>* arg) { return std::atan(arg[0]); }
	std::complex<double> _sinhC(std::complex<double>* arg) { return std::sinh(arg[0]); }
	std::complex<double> _coshC(std::complex<double>* arg) { return std::cosh(arg[0]); }
	std::complex<double> _tanhC(std::complex<double>* arg) { return std::tanh(arg[0]); }

	// helpers for inline.
	std::complex<double> dummyCC(std::complex<double> * args) { return std::complex<double>(0, 0); }
	double dummyRR(double args) { return 0; }



	// helpers, no derived object
	double logRR(double x) { return std::log(x); }
	std::complex<double> logCC(std::complex<double> *  x) { return std::log(x[0]); }
	double log2RR(double x) { return std::log2(x); }
	std::complex<double> log2CC(std::complex<double> *  x) { return std::log(x[0]) / log(2); }
	double log10RR(double x) { return std::log10(x); }
	std::complex<double> log10CC(std::complex<double> *  x) { return std::log10(x[0]); }
	double expRR(double x) { return std::exp(x); }
	std::complex<double> expCC(std::complex<double> *  x) { return std::exp(x[0]); }
	double powRR(double x, double y) { return std::pow(x, y); }
	std::complex<double> powCC(std::complex<double> *  x) { return std::pow(x[0], x[1]); }
	double atan2RR(double x, double y) { return std::atan2(x, y); }
	double hypotRR(double x, double y) { return std::hypot(x, y); }
	std::complex<double> sqrtRC(double  * x) { return std::sqrt(std::complex<double>(x[0], 0)); }
	std::complex<double> sqrtCC(std::complex<double> *  x) { return std::sqrt(x[0]); }


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
		ctx->addObject(":", std::make_shared<MpOperationTernaryHelper>());
		ctx->addObject("=", std::make_shared<MpOperationAssignment>());
		ctx->addObject("isfinite", std::make_shared<MpOperationIsFinite>());
		ctx->addObject("isinf", std::make_shared<MpOperationIsInfinite>());
		ctx->addObject("isnan", std::make_shared<MpOperationIsNan>());
		ctx->addObject("real", std::make_shared<MpOperationGetReal>());
		ctx->addObject("imag", std::make_shared<MpOperationGetImag>());
		ctx->addObject("min", std::make_shared<MpOperationMin>());
		ctx->addObject("max", std::make_shared<MpOperationMax>());
		ctx->addObject("=", std::make_shared<MpOperationAssignment>());
		ctx->addObject("%", std::make_shared<MpOperationModulo>());
		ctx->addObject("-", std::make_shared<MpOperationNeg>());
		ctx->addObject("!", std::make_shared<MpOperationNot>());
		ctx->addObject("sqrt", std::make_shared<MpOperationSqrt>());
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
		ctx->addObject("sin", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::sin));
		ctx->addObject("cos", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::cos));
		ctx->addObject("tan", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::tan));
		ctx->addObject("sinh", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::sinh));
		ctx->addObject("cosh", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::cosh));
		ctx->addObject("tanh", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::tanh));
		ctx->addObject("asin", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::asin));
		ctx->addObject("acos", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::acos));
		ctx->addObject("atan", std::make_shared<MpOperationTrigonometrie>(MpOperationTrigonometrie::atan));
		ctx->addObject("sqrtc", std::make_shared<MpOperationFunc>(1, MpOperationFlags::OpFlagDReturnsC, VPTR(sqrtRC), VPTR(sqrtCC)));
		ctx->addObject("log", std::make_shared<MpOperationFunc>(1, MpOperationFlags::OpFlagNone, VPTR(logRR), VPTR(logCC)));
		ctx->addObject("log2", std::make_shared<MpOperationFunc>(1, MpOperationFlags::OpFlagNone, VPTR(log2RR), VPTR(log2CC)));
		ctx->addObject("log10", std::make_shared<MpOperationFunc>(1, MpOperationFlags::OpFlagNone, VPTR(log10RR), VPTR(log10CC)));
		ctx->addObject("exp", std::make_shared<MpOperationFunc>(1, MpOperationFlags::OpFlagNone, VPTR(expRR), VPTR(expCC)));
		ctx->addObject("pow", std::make_shared<MpOperationFunc>(2, MpOperationFlags::OpFlagNone, VPTR(powRR), VPTR(powCC)));
		ctx->addObject("atan2", std::make_shared<MpOperationFunc>(2, MpOperationFlags::OpFlagNone, VPTR(atan2RR), nullptr));
		ctx->addObject("hypot", std::make_shared<MpOperationFunc>(2, MpOperationFlags::OpFlagNone, VPTR(hypotRR), nullptr));
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
			if (!fnC_)
				return opt->_errorReporter->onError(kErrorInvalidArgument, node->getPosition(),
					"No complex function available.");

			node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
			b_returns_complex = !hasFlag(MpOperationFlags::OpFlagCReturnsD);
		}
		else
		{
			if (fnD_)
			{
				node->removeNodeFlags(AstNodeFlags::kAstTakesComplex);
				b_returns_complex = hasFlag(MpOperationFlags::OpFlagDReturnsC);
			}
			else if (fnC_)
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
		if (b_all_imm && !hasFlag(MpOperationFlags::OpFlagHasState))
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
		}
		else
		{
			fnD_ = fn;
		}
	}

	double MpOperationFunc::evaluateDRetD(double * args)
	{
		if (!fnD_)
		{
			throw std::runtime_error("Function does not exist.");
		}

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
	double isfiniteRR(double args) { return std::isfinite(args) ? 1.0 : 0.0; }
	std::complex<double> isfiniteCC(std::complex<double>* args)
	{
		return std::complex<double>(std::isfinite(args[0].real()) ? 1.0 : 0.0, std::isfinite(args[0].imag()) ? 1.0 : 0.0);
	}

	MpOperationIsFinite::MpOperationIsFinite() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(isfiniteRR), VPTR(isfiniteCC), nullptr, nullptr)
	{}
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

	// MpOperationIsInFinite
	double isinfRR(double args) { return std::isinf(args) ? 1.0 : 0.0; }
	std::complex<double> isinfCC(std::complex<double>* args)
	{
		return std::complex<double>(std::isinf(args[0].real()) ? 1.0 : 0.0, std::isinf(args[0].imag()) ? 1.0 : 0.0);
	}

	MpOperationIsInfinite::MpOperationIsInfinite() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(isinfRR), VPTR(isinfCC), nullptr, nullptr)
	{}

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

	// MpOperationIsNan
	double isnanRR(double args) { return std::isnan(args) ? 1.0 : 0.0; }
	std::complex<double> isnanCC(std::complex<double>* args)
	{
		return std::complex<double>(std::isnan(args[0].real()) ? 1.0 : 0.0, std::isnan(args[0].imag()) ? 1.0 : 0.0);
	}

	MpOperationIsNan::MpOperationIsNan() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(isnanRR), VPTR(isnanCC), nullptr, nullptr)
	{}

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

	// MpOperationGetReal
	double realCR(std::complex<double>* args) { return args->real(); }

	MpOperationGetReal::MpOperationGetReal() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, nullptr, VPTR(realCR), nullptr, nullptr)
	{}

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


	// MpOperationGetImag
	double imagCR(std::complex<double>* args) { return args->imag(); }

	MpOperationGetImag::MpOperationGetImag() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, nullptr, VPTR(imagCR), nullptr, nullptr)
	{}

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

	// Square root
	double sqrtRR(double args) { return std::sqrt(args); }

	MpOperationSqrt::MpOperationSqrt() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(sqrtRR), nullptr, nullptr, nullptr)
	{}

	JitVar MpOperationSqrt::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var = jc->onNode(node->getAt(0));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAGS::FLAG_NONE);
		if (var.isXmm())
			jc->cc->sqrtsd(result.getXmm(), var.getXmm());
		else
			jc->cc->sqrtsd(result.getXmm(), var.getMem());
		return result;
	}

	// Negation
	double negRR(double args) { return -args; }
	std::complex<double> negCC(std::complex<double>* args) { return -args[0]; }

	MpOperationNeg::MpOperationNeg() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(negRR), VPTR(negCC), nullptr, nullptr)
	{
		priority_ = 3;
	}

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
	double notRR(double args) { return args == 0 ? 1.0 : 0.0; }
	std::complex<double> notCC(std::complex<double>* args) { return std::complex<double>(args[0] == std::complex<double>(0, 0) ? 1.0 : 0.0, 0.0); }

	MpOperationNot::MpOperationNot() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(notRR), VPTR(notCC), nullptr, nullptr)
	{
		priority_ = 3;
	}

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

	// Conjugate
	std::complex<double> conjugCC(std::complex<double>* args) { return std::complex<double>(args->real(), -args->imag()); }

	MpOperationConjug::MpOperationConjug() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, VPTR(conjugCC), nullptr, nullptr)
	{}

	JitVar MpOperationConjug::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar tmp = jc->onNode(node->getAt(0));
		JitVar result = jc->registerVarComplex(tmp, !node->getAt(0)->returnsComplex());
		jc->cc->pxor(result.getXmm(), jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000)).getMem());
		return result;
	}

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
	std::complex<double> recipCC(std::complex<double>* args) { return 1.0 / args[0]; }
	double recipRR(double args) { return 1.0 / args; }

	MpOperationRecip::MpOperationRecip() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(recipRR), VPTR(recipCC), nullptr, nullptr)
	{}

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


	// MpOperationTrigonometrie

	MpOperationTrigonometrie::MpOperationTrigonometrie(uint32_t type) :
		MpOperationFuncAsm(1, OpFlagHasAsm, VPTR(dummyRR), VPTR(dummyCC), nullptr, nullptr),
		type_(type)
	{}

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
		case trigonometrieFunc::sin: return _sin(args[0]);
		case trigonometrieFunc::cos: return _cos(args[0]);
		case trigonometrieFunc::tan: return _tan(args[0]);
		case trigonometrieFunc::asin: return _asin(args[0]);
		case trigonometrieFunc::acos: return _acos(args[0]);
		case trigonometrieFunc::atan: return _atan(args[0]);
		case trigonometrieFunc::sinh: return _sinh(args[0]);
		case trigonometrieFunc::cosh: return _cosh(args[0]);
		case trigonometrieFunc::tanh: return _tanh(args[0]);
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
	double signbitRR(double args) { return std::signbit(args) ? 1.0 : 0.0; }

	MpOperationSignBit::MpOperationSignBit() :
		MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(signbitRR), nullptr, nullptr, nullptr)
	{}
	JitVar MpOperationSignBit::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar result(jc->cc->newXmmSd(), JitVar::FLAG_NONE);
		jc->cc->pshufd(result.getXmm(), jc->registerVar(var).getXmm(), asmjit::x86::shufImm(3, 2, 1, 1));
		jc->cc->psrad(result.getXmm(), 31);
		jc->cc->andpd(result.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return result;
	}

	// Copy sign
	double copysignRR(double args0, double args1) { return std::copysign(args0, args1); }

	MpOperationCopySign::MpOperationCopySign() : MpOperationFuncAsm(2, MpOperationFlags::OpFlagNone, VPTR(copysignRR), nullptr, nullptr, nullptr)
	{}

	JitVar MpOperationCopySign::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar vl = jc->writableVar(jc->onNode(node->getAt(0)));
		JitVar vr = jc->writableVar(jc->onNode(node->getAt(1)));

		jc->cc->andpd(vl.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x7FFFFFFFFFFFFFFF)).getMem());
		jc->cc->andpd(vr.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->orpd(vl.getXmm(), vr.getXmm());

		return vl;
	}


	// Average
	double avgRR(double args0, double args1) { return (args0 + args1) * 0.5; }
	std::complex<double> avgCC(std::complex<double>* args) { return (args[0] + args[1]) * 0.5; }

	MpOperationAvg::MpOperationAvg() : MpOperationFuncAsm(2, MpOperationFlags::OpFlagNone, VPTR(avgRR), VPTR(avgCC), nullptr, nullptr)
	{}

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
	
	// Absolute
	double absRR(double args) { return std::abs(args); }
	double absCR(std::complex<double>* args) { return std::abs(args[0]); }

	MpOperationAbs::MpOperationAbs() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, VPTR(absRR), VPTR(absCR), nullptr, nullptr)
	{}

	JitVar MpOperationAbs::compile(JitCompiler * jc, AstNode * node)
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar result;
		if (node->takesComplex())
		{
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



	// round
	double roundRR(double args)
	{
		double y = ::floor(args);
		return y + (args - y >= 0.5 ? double(1.0) : double(0.0));
	}

	MpOperationRound::MpOperationRound() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(roundRR), nullptr, nullptr, nullptr)
	{}

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
	
	// roundeven
	double roundevenRR(double args) { return std::rint(args); }

	MpOperationRoundEven::MpOperationRoundEven() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(roundevenRR), nullptr, nullptr, nullptr)
	{}

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



	// trunc
	double truncRR(double args) { return std::trunc(args); }

	MpOperationTrunc::MpOperationTrunc() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(truncRR), nullptr, nullptr, nullptr)
	{}

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


	// frac
	double fracRR(double args) { return args - std::floor(args); }

	MpOperationFrac::MpOperationFrac() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(fracRR), nullptr, nullptr, nullptr)
	{}

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


	// floor
	double floorRR(double args) { return std::floor(args); }

	MpOperationFloor::MpOperationFloor() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(floorRR), nullptr, nullptr, nullptr)
	{}

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


	// ceil
	double ceilRR(double args) { return std::ceil(args); }

	MpOperationcCeil::MpOperationcCeil() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, VPTR(ceilRR), nullptr, nullptr, nullptr)
	{}

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

		if (lIsImm && rIsImm && !hasFlag(MpOperationFlags::OpFlagHasState))
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
		else if (lIsImm && !hasFlag(MpOperationFlags::OpFlagHasState))
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
		else if (rIsImm && !hasFlag(MpOperationFlags::OpFlagHasState))
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
		ternaryNode->mpOp_ = opt->_symbols->at(std::make_pair("?", 2)).get();

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

	JitVar MpOperationTernaryHelper::compile(JitCompiler * jc, AstNode * node)
	{
		throw std::runtime_error("MpOperationHelper: should never be reached");
	}

	uint32_t MpOperationTernaryHelper::optimize(AstOptimizer *opt, AstNode *node)
	{
		return kErrorInvalidState;
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
