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
#include <mathpresso/mpeval_p.h>
#include <mathpresso/mpoptimizer_p.h>
#include <asmjit/x86/x86operand.h>
#include <asmjit/x86/x86inst.h>

#include <complex>

namespace mathpresso
{

#define VPTR(function) reinterpret_cast<void*>(function)

#ifdef _REALREWORK
	double sinRR(double * arg) { return std::sin(arg[0]); }
	double cosRR(double * arg) { return std::cos(arg[0]); }
	double tanRR(double * arg) { return std::tan(arg[0]); }
	double asinRR(double * arg) { return std::asin(arg[0]); }
	double acosRR(double * arg) { return std::acos(arg[0]); }
	double atanRR(double * arg) { return std::atan(arg[0]); }
	double sinhRR(double * arg) { return std::sinh(arg[0]); }
	double coshRR(double * arg) { return std::cosh(arg[0]); }
	double tanhRR(double * arg) { return std::tanh(arg[0]); }

	double logRR(double* x) { return std::log(x[0]); }
	double log2RR(double * x) { return std::log2(x[0]); }
	double log10RR(double * x) { return std::log10(x[0]); }
	double powRR(double * x) { return std::pow(x[0], x[1]); }
	double expRR(double * x) { return std::exp(x[0]); }
	double atan2RR(double * x) { return std::atan2(x[0], x[1]); }
	double hypotRR(double * x) { return std::hypot(x[0], x[1]); }

	double isfiniteRR(double * args) { return std::isfinite(args[0]) ? 1.0 : 0.0; }
	double isinfRR(double * args) { return std::isinf(args[0]) ? 1.0 : 0.0; }
	double isnanRR(double *  args) { return std::isnan(args[0]) ? 1.0 : 0.0; }
	double sqrtRR(double * args) { return std::sqrt(args[0]); }
	double negRR(double * args) { return -args[0]; }
	double notRR(double * args) { return args[0] == 0 ? 1.0 : 0.0; }
	double recipRR(double * args) { return 1.0 / args[0]; }
	double signbitRR(double * args) { return std::signbit(args[0]) ? 1.0 : 0.0; }
	double copysignRR(double * args) { return std::copysign(args[0], args[1]); }
	double avgRR(double * args) { return (args[0] + args[1]) * 0.5; }
	double absRR(double * args) { return std::abs(args[0]); }
	double roundRR(double * args) { return std::floor(args[0] + .5); }
	double roundevenRR(double * args) { return std::rint(args[0]); }
	double truncRR(double * args) { return std::trunc(args[0]); }
	double fracRR(double * args) { return args[0] - std::floor(args[0]); }
	double floorRR(double * args) { return std::floor(args[0]); }
	double ceilRR(double * args) { return std::ceil(args[0]); }
#else
	double sinRR(double arg) { return std::sin(arg); }
	double cosRR(double arg) { return std::cos(arg); }
	double tanRR(double arg) { return std::tan(arg); }
	double asinRR(double arg) { return std::asin(arg); }
	double acosRR(double arg) { return std::acos(arg); }
	double atanRR(double arg) { return std::atan(arg); }
	double sinhRR(double arg) { return std::sinh(arg); }
	double coshRR(double arg) { return std::cosh(arg); }
	double tanhRR(double arg) { return std::tanh(arg); }

	double logRR(double x) { return std::log(x); }
	double log2RR(double x) { return std::log2(x); }
	double log10RR(double x) { return std::log10(x); }
	double expRR(double x) { return std::exp(x); }
	double powRR(double x, double y) { return std::pow(x, y); }
	double atan2RR(double x, double y) { return std::atan2(x, y); }
	double hypotRR(double x, double y) { return std::hypot(x, y); }
	double isfiniteRR(double args) { return std::isfinite(args) ? 1.0 : 0.0; }
	double isinfRR(double args) { return std::isinf(args) ? 1.0 : 0.0; }
	double isnanRR(double args) { return std::isnan(args) ? 1.0 : 0.0; }
	double sqrtRR(double args) { return std::sqrt(args); }
	double negRR(double args) { return -args; }
	double notRR(double args) { return args == 0 ? 1.0 : 0.0; }
	double recipRR(double args) { return 1.0 / args; }
	double signbitRR(double args) { return std::signbit(args) ? 1.0 : 0.0; }
	double copysignRR(double args0, double args1) { return std::copysign(args0, args1); }
	double avgRR(double args0, double args1) { return (args0 + args1) * 0.5; }
	double absRR(double args) { return std::abs(args); }
	double roundRR(double args) { return std::floor(args + .5); }
	double roundevenRR(double args) { return std::rint(args); }
	double truncRR(double args) { return std::trunc(args); }
	double fracRR(double args) { return args - std::floor(args); }
	double floorRR(double args) { return std::floor(args); }
	double ceilRR(double args) { return std::ceil(args); }
#endif

	// helpers, no derived object
	std::complex<double> sinCC(std::complex<double>* arg) { return std::sin(arg[0]); }
	std::complex<double> cosCC(std::complex<double>* arg) { return std::cos(arg[0]); }
	std::complex<double> tanCC(std::complex<double>* arg) { return std::tan(arg[0]); }
	std::complex<double> asinCC(std::complex<double>* arg) { return std::asin(arg[0]); }
	std::complex<double> acosCC(std::complex<double>* arg) { return std::acos(arg[0]); }
	std::complex<double> atanCC(std::complex<double>* arg) { return std::atan(arg[0]); }
	std::complex<double> sinhCC(std::complex<double>* arg) { return std::sinh(arg[0]); }
	std::complex<double> coshCC(std::complex<double>* arg) { return std::cosh(arg[0]); }
	std::complex<double> tanhCC(std::complex<double>* arg) { return std::tanh(arg[0]); }
	std::complex<double> logCC(std::complex<double> *  x) { return std::log(x[0]); }
	std::complex<double> log2CC(std::complex<double> *  x) { return std::log(x[0]) / log(2); }
	std::complex<double> log10CC(std::complex<double> *  x) { return std::log10(x[0]); }
	std::complex<double> expCC(std::complex<double> *  x) { return std::exp(x[0]); }
	std::complex<double> powCC(std::complex<double> *  x) { return std::pow(x[0], x[1]); }
	std::complex<double> sqrtRC(double  * x) { return std::sqrt(std::complex<double>(x[0], 0)); }
	std::complex<double> sqrtCC(std::complex<double> *  x) { return std::sqrt(x[0]); }
	double absCR(std::complex<double>* args) { return std::abs(args[0]); }


	uint32_t addBuiltinMpObjects(Context * ctx)
	{
		ctx->addObject("+", std::make_shared<MpOperationAdd<double>>());
		ctx->addObject("+", std::make_shared<MpOperationAdd<std::complex<double>>>());
		ctx->addObject("-", std::make_shared<MpOperationSub<double>>());
		ctx->addObject("-", std::make_shared<MpOperationSub<std::complex<double>>>());
		ctx->addObject("*", std::make_shared<MpOperationMul<double>>());
		ctx->addObject("*", std::make_shared<MpOperationMul<std::complex<double>>>());
		ctx->addObject("/", std::make_shared<MpOperationDiv<double>>());
		ctx->addObject("/", std::make_shared<MpOperationDiv<std::complex<double>>>());
		ctx->addObject("==", std::make_shared<MpOperationEq<double>>());
		ctx->addObject("==", std::make_shared<MpOperationEq<std::complex<double>>>());
		ctx->addObject("!=", std::make_shared<MpOperationNe<double>>());
		ctx->addObject("!=", std::make_shared<MpOperationNe<std::complex<double>>>());
		ctx->addObject(">=", std::make_shared<MpOperationGe>());
		ctx->addObject(">", std::make_shared<MpOperationGt>());
		ctx->addObject("<=", std::make_shared<MpOperationLe>());
		ctx->addObject("<", std::make_shared<MpOperationLt>());
		ctx->addObject("_ternary_", std::make_shared<MpOperationTernary<double>>());
		ctx->addObject("_ternary_", std::make_shared<MpOperationTernary<std::complex<double>>>());
		ctx->addObject("=", std::make_shared<MpOperationAssignment<double>>());
		ctx->addObject("=", std::make_shared<MpOperationAssignment<std::complex<double>>>());
		ctx->addObject("isfinite", std::make_shared<MpOperationIsFinite<double>>());
		ctx->addObject("isfinite", std::make_shared<MpOperationIsFinite<std::complex<double>>>());
		ctx->addObject("isinf", std::make_shared<MpOperationIsInfinite<double>>());
		ctx->addObject("isinf", std::make_shared<MpOperationIsInfinite<std::complex<double>>>());
		ctx->addObject("isnan", std::make_shared<MpOperationIsNan<double>>());
		ctx->addObject("isnan", std::make_shared<MpOperationIsNan<std::complex<double>>>());

		ctx->addObject("-", std::make_shared<MpOperationNeg<double>>());
		ctx->addObject("-", std::make_shared<MpOperationNeg<std::complex<double>>>());
		ctx->addObject("abs", std::make_shared<MpOperationAbs>());
		ctx->addObject("abs", std::make_shared<MpOperationFunc<double, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(absCR)));
		ctx->addObject("avg", std::make_shared<MpOperationAvg<double>>());
		ctx->addObject("avg", std::make_shared<MpOperationAvg<std::complex<double>>>());
		ctx->addObject("recip", std::make_shared<MpOperationRecip<double>>());
		ctx->addObject("recip", std::make_shared<MpOperationRecip<std::complex<double>>>());
		ctx->addObject("!", std::make_shared<MpOperationNot<double>>());
		ctx->addObject("!", std::make_shared<MpOperationNot<std::complex<double>>>());

		ctx->addObject("real", std::make_shared<MpOperationGetReal>());
		ctx->addObject("imag", std::make_shared<MpOperationGetImag>());
		ctx->addObject("min", std::make_shared<MpOperationMin>());
		ctx->addObject("max", std::make_shared<MpOperationMax>());
		ctx->addObject("%", std::make_shared<MpOperationModulo>());
		ctx->addObject("sqrt", std::make_shared<MpOperationSqrt>());
		ctx->addObject("conjug", std::make_shared<MpOperationConjug>());
		ctx->addObject("signbit", std::make_shared<MpOperationSignBit>());
		ctx->addObject("copysign", std::make_shared<MpOperationCopySign>());
		ctx->addObject("round", std::make_shared<MpOperationRound>());
		ctx->addObject("roundeven", std::make_shared<MpOperationRoundEven>());
		ctx->addObject("floor", std::make_shared<MpOperationFloor>());
		ctx->addObject("ceil", std::make_shared<MpOperationcCeil>());
		ctx->addObject("frac", std::make_shared<MpOperationFrac>());
		ctx->addObject("trunc", std::make_shared<MpOperationTrunc>());

		ctx->addObject("sin", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(sinRR)));
		ctx->addObject("sin", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(sinCC)));
		ctx->addObject("cos", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(cosRR)));
		ctx->addObject("cos", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(cosCC)));
		ctx->addObject("tan", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(tanRR)));
		ctx->addObject("tan", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(tanCC)));
		ctx->addObject("sinh", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(sinhRR)));
		ctx->addObject("sinh", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(sinhCC)));
		ctx->addObject("cosh", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(coshRR)));
		ctx->addObject("cosh", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(coshCC)));
		ctx->addObject("tanh", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(tanhRR)));
		ctx->addObject("tanh", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(tanhCC)));
		ctx->addObject("asin", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(asinRR)));
		ctx->addObject("asin", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(asinCC)));
		ctx->addObject("acos", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(acosRR)));
		ctx->addObject("acos", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(acosCC)));
		ctx->addObject("atan", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(atanRR)));
		ctx->addObject("atan", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(atanCC)));
		ctx->addObject("sqrtc", std::make_shared<MpOperationFunc<std::complex<double>, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(sqrtRC)));
		ctx->addObject("sqrtc", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(sqrtCC)));
		ctx->addObject("log", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(logRR)));
		ctx->addObject("log", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(logCC)));
		ctx->addObject("log2", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(log2RR)));
		ctx->addObject("log2", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(log2CC)));
		ctx->addObject("log10", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(log10RR)));
		ctx->addObject("log10", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(log10CC)));
		ctx->addObject("exp", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 1, VPTR(expRR)));
		ctx->addObject("exp", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 1, VPTR(expCC)));
		ctx->addObject("pow", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 2, VPTR(powRR)));
		ctx->addObject("pow", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 2, VPTR(powCC)));
		ctx->addObject("atan2", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 2, VPTR(atan2RR)));
		ctx->addObject("hypot", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 2, VPTR(hypotRR)));

		ctx->addObject("_none_", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpFlagNone, 0, nullptr));
		ctx->addObject("_none_", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpFlagNone, 0, nullptr));

		ctx->addObject("?", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpIsRighttoLeft, 2, nullptr, 15));
		ctx->addObject("?", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpIsRighttoLeft, 2, nullptr, 15));
		ctx->addObject(":", std::make_shared<MpOperationFunc<double, double>>(MpOperationFlags::OpIsRighttoLeft, 2, nullptr, 15));
		ctx->addObject(":", std::make_shared<MpOperationFunc<std::complex<double>, std::complex<double>>>(MpOperationFlags::OpIsRighttoLeft, 2, nullptr, 15));

		ctx->addConstant("NaN", mpGetNan());
		ctx->addConstant("INF", mpGetInf());
		ctx->addConstant("PI", 3.14159265358979323846);
		ctx->addConstant("E", 2.7182818284590452354);
		ctx->addConstant("i", { 0, 1 });
		return 0;
	}

	void Signature::init(type retType, std::vector<param> params, uint32_t flags)
	{
		return_type_ = retType;
		parameters_ = params;
		flags_ = flags;
	}

	template<>
	MpOperationFunc<double, double>::MpOperationFunc(uint32_t flags, size_t numargs, void * fnPtr, uint32_t priority) :
		MpOperation(Signature(numargs, Signature::type::real, flags), priority),
		fnPtr_(fnPtr)
	{
	}

	template<>
	MpOperationFunc<std::complex<double>, double>::MpOperationFunc(uint32_t flags, size_t numargs, void * fnPtr, uint32_t priority) :
		MpOperation(Signature(Signature::type::complex, { numargs,{ Signature::type::real } }, flags), priority),
		fnPtr_(fnPtr)
	{
	}

	template<>
	MpOperationFunc<double, std::complex<double>>::MpOperationFunc(uint32_t flags, size_t numargs, void * fnPtr, uint32_t priority) :
		MpOperation(Signature(Signature::type::real, { numargs, {Signature::type::complex} }, flags), priority),
		fnPtr_(fnPtr)
	{
	}

	template<>
	MpOperationFunc<std::complex<double>, std::complex<double>>::MpOperationFunc(uint32_t flags, size_t numargs, void * fnPtr, uint32_t priority) :
		MpOperation(Signature(numargs, Signature::type::complex, flags), priority),
		fnPtr_(fnPtr)
	{
	}


	template<>
	inline JitVar MpOperationFunc<double, double>::compile(JitCompiler * jc, AstNode * node) const
	{
		asmjit::X86Xmm result = jc->cc->newXmmSd();
		asmjit::X86Xmm args[8];
		if (!fnPtr_)
		{
			// Should never happen, as the optimizer should have taken care of that. Remove later
			throw std::runtime_error("Implementation error!");
		}

		for (size_t i = 0; i < nargs_; i++)
		{
			args[i] = jc->registerVar(jc->onNode(node->getAt(i))).getXmm();
		}
		jc->inlineCall<double, double>(result, args, nargs_, fnPtr_);

		return JitVar(result, false);
	}
	template<>
	inline JitVar MpOperationFunc<double, std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		asmjit::X86Xmm result = jc->cc->newXmmSd();
		asmjit::X86Xmm args[8];
		if (!fnPtr_)
		{
			// Should never happen, as the optimizer should have taken care of that. Remove later
			throw std::runtime_error("Implementation error!");
		}

		for (size_t i = 0; i < nargs_; i++)
		{
			args[i] = jc->registerVarComplex(jc->onNode(node->getAt(i)), !node->getAt(i)->returnsComplex()).getXmm();
		}

		jc->inlineCall<double, std::complex<double>>(result, args, nargs_, fnPtr_);

		return JitVar(result, false);
	}
	template<>
	inline JitVar MpOperationFunc<std::complex<double>, double>::compile(JitCompiler * jc, AstNode * node) const
	{
		asmjit::X86Xmm result = jc->cc->newXmmPd();
		asmjit::X86Xmm args[8];
		if (!fnPtr_)
		{
			// Should never happen, as the optimizer should have taken care of that. Remove later
			throw std::runtime_error("Implementation error!");
		}

		for (size_t i = 0; i < nargs_; i++)
		{
			args[i] = jc->registerVar(jc->onNode(node->getAt(i))).getXmm();
		}
		jc->inlineCall<std::complex<double>, double>(result, args, nargs_, fnPtr_);

		return JitVar(result, false);
	}
	template<>
	inline JitVar MpOperationFunc<std::complex<double>, std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		asmjit::X86Xmm result = jc->cc->newXmmPd();
		asmjit::X86Xmm args[8];
		if (!fnPtr_)
		{
			// Should never happen, as the optimizer should have taken care of that. Remove later
			throw std::runtime_error("Implementation error!");
		}

		for (size_t i = 0; i < nargs_; i++)
		{
			args[i] = jc->registerVarComplex(jc->onNode(node->getAt(i)), !node->getAt(i)->returnsComplex()).getXmm();
		}

		jc->inlineCall<std::complex<double>, std::complex<double>>(result, args, nargs_, fnPtr_);

		return JitVar(result, false);
	}

	template<typename RET, typename PARAM>
	uint32_t MpOperationFunc<RET, PARAM>::optimize(AstOptimizer * opt, AstNode * node) const
	{
		bool b_all_imm = true;

		// Gather Information about the child-nodes.
		for (size_t i = 0; i < node->getLength(); i++)
		{
			b_all_imm &= node->getAt(i)->isImm();
		}

		// optimize all-immediate calls:
		if (b_all_imm && !hasFlag(MpOperationFlags::OpFlagHasState))
		{
			AstImm* ret = opt->getAst()->newNode<AstImm>(0);

			std::vector<PARAM> args;
			for (size_t i = 0; i < nargs_; i++)
			{
				args.push_back((static_cast<AstImm*>(node->getAt(i)))->getValue<PARAM>());
			}

			ret->setValue(evaluate(args.data()));

			node->getParent()->replaceNode(node, ret);
			opt->getAst()->deleteNode(node);
			node = ret;
		}
		return ErrorCode::kErrorOk;
	}

	template<>
	double MpOperationFunc<double, double>::evaluate(double * args) const
	{
		if (!fnPtr_)
		{
			throw std::runtime_error("Function does not exist.");
		}
#ifdef _REALREWORK
		return ((mpFuncDtoD)fnPtr_)(args);
#else
		switch (nargs_)
		{
			case 0: return ((Arg0Func)fnPtr_)();
			case 1: return ((Arg1Func)fnPtr_)(args[0]);
			case 2: return ((Arg2Func)fnPtr_)(args[0], args[1]);
			case 3: return ((Arg3Func)fnPtr_)(args[0], args[1], args[2]);
			case 4: return ((Arg4Func)fnPtr_)(args[0], args[1], args[2], args[3]);
			case 5: return ((Arg5Func)fnPtr_)(args[0], args[1], args[2], args[3], args[4]);
			case 6: return ((Arg6Func)fnPtr_)(args[0], args[1], args[2], args[3], args[4], args[5]);
			case 7: return ((Arg7Func)fnPtr_)(args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
			case 8: return ((Arg8Func)fnPtr_)(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]);
			default:
				throw std::runtime_error("Too many arguments.");
		}
#endif // _REALREWORK
	}
	template<>
	std::complex<double> MpOperationFunc<std::complex<double>, double>::evaluate(double * args) const
	{
		if (!fnPtr_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpDtoC)fnPtr_)(args);
	}
	template<>
	double MpOperationFunc<double, std::complex<double>>::evaluate(std::complex<double>* args) const
	{
		if (!fnPtr_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpCtoD)fnPtr_)(args);
	}
	template<>
	std::complex<double> MpOperationFunc<std::complex<double>, std::complex<double>>::evaluate(std::complex<double>* args) const
	{
		if (!fnPtr_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpCtoC)fnPtr_)(args);
	}

	// MpOperationIsFinite
	std::complex<double> isfiniteCC(std::complex<double>* args)
	{
		return std::complex<double>(std::isfinite(args[0].real()) ? 1.0 : 0.0, std::isfinite(args[0].imag()) ? 1.0 : 0.0);
	}
	template<>
	MpOperationIsFinite<double>::MpOperationIsFinite() :
		MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(isfiniteRR))
	{
	}
	template<>
	MpOperationIsFinite<std::complex<double>>::MpOperationIsFinite() :
		MpOperationFunc<std::complex<double>, std::complex<double>>(MpOperationFlags::OpFlagNone, 1, VPTR(isfiniteCC))
	{
	}

	template<>
	JitVar MpOperationIsFinite<double>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVar(var);
		jc->cc->orpd(var.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->cmpsd(var.getXmm(), jc->getConstantU64(0).getMem(), int(asmjit::x86::kCmpLE));
		jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return var;
	}

	template<>
	JitVar MpOperationIsFinite<std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->orpd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0x8000000000000000), MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->cmppd(var.getXmm(), jc->getConstantD64(std::complex<double>(0.0, 0.0)).getMem(), int(asmjit::x86::kCmpLE));
		jc->cc->andpd(var.getXmm(), jc->getConstantD64(std::complex<double>(1.0, 1.0)).getMem());
		return var;
	}


	// MpOperationIsInFinite
	std::complex<double> isinfCC(std::complex<double>* args)
	{
		return std::complex<double>(std::isinf(args[0].real()) ? 1.0 : 0.0, std::isinf(args[0].imag()) ? 1.0 : 0.0);
	}
	template<>
	MpOperationIsInfinite<double>::MpOperationIsInfinite() :
		MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(isinfRR))
	{
	}
	template<>
	MpOperationIsInfinite<std::complex<double>>::MpOperationIsInfinite() :
		MpOperationFunc<std::complex<double>, std::complex<double>>(MpOperationFlags::OpFlagNone, 1, VPTR(isinfCC))
	{
	}

	template<>
	JitVar MpOperationIsInfinite<double>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVar(var);
		jc->cc->orpd(var.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->cmpsd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0xFFF0000000000000)).getMem(), int(asmjit::x86::kCmpEQ));
		jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return var;
	}

	template<>
	JitVar MpOperationIsInfinite<std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->orpd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0x8000000000000000), MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->cmppd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0xFFF0000000000000), MATHPRESSO_UINT64_C(0xFFF0000000000000)).getMem(), int(asmjit::x86::kCmpEQ));
		jc->cc->andpd(var.getXmm(), jc->getConstantD64(std::complex<double>(1.0, 1.0)).getMem());
		return var;
	}


	// MpOperationIsNan	
	std::complex<double> isnanCC(std::complex<double>* args)
	{
		return std::complex<double>(std::isnan(args[0].real()) ? 1.0 : 0.0, std::isnan(args[0].imag()) ? 1.0 : 0.0);
	}
	template<>
	MpOperationIsNan<double>::MpOperationIsNan() :
		MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(isnanRR))
	{
	}
	template<>
	MpOperationIsNan<std::complex<double>>::MpOperationIsNan() :
		MpOperationFunc<std::complex<double>, std::complex<double>>(MpOperationFlags::OpFlagNone, 1, VPTR(isnanCC))
	{
	}

	template<>
	JitVar MpOperationIsNan<double>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVar(var);
		jc->cc->cmpsd(var.getXmm(), var.getXmm(), int(asmjit::x86::kCmpEQ));
		jc->cc->andnpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return var;
	}

	template<>
	JitVar MpOperationIsNan<std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->cmppd(var.getXmm(), var.getXmm(), int(asmjit::x86::kCmpEQ));
		jc->cc->andnpd(var.getXmm(), jc->getConstantD64(std::complex<double>(1.0, 1.0)).getMem());
		return var;
	}

	// MpOperationGetReal
	double realCR(std::complex<double>* args)
	{
		return args->real();
	}

	MpOperationGetReal::MpOperationGetReal() :
		MpOperationFunc<double, std::complex<double>>(MpOperationFlags::OpFlagNone, 1, VPTR(realCR))
	{
	}

	JitVar MpOperationGetReal::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar varRet(jc->cc->newXmmSd(), false);;
		jc->cc->xorpd(varRet.getXmm(), varRet.getXmm());
		if (var.isXmm())
		{
			jc->cc->movsd(varRet.getXmm(), var.getXmm());
		}
		else
		{
			jc->cc->movsd(varRet.getXmm(), var.getMem());
		}
		return varRet;
	}


	// MpOperationGetImag
	double imagCR(std::complex<double>* args)
	{
		return args->imag();
	}

	MpOperationGetImag::MpOperationGetImag() :
		MpOperationFunc<double, std::complex<double>>(MpOperationFlags::OpFlagNone, 1, VPTR(imagCR))
	{
	}

	JitVar MpOperationGetImag::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar varRet(jc->cc->newXmmSd(), false);;
		var = jc->registerVarComplex(var, !node->getAt(0)->hasNodeFlag(AstNodeFlags::kAstReturnsComplex));
		jc->cc->xorpd(varRet.getXmm(), varRet.getXmm());
		jc->cc->shufpd(var.getXmm(), varRet.getXmm(), asmjit::x86::shufImm(0, 1));
		return var;
	}

	// Square root
	MpOperationSqrt::MpOperationSqrt() :
		MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(sqrtRR))
	{
	}

	JitVar MpOperationSqrt::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		JitVar result(jc->cc->newXmmSd(), false);
		if (var.isXmm())
			jc->cc->sqrtsd(result.getXmm(), var.getXmm());
		else
			jc->cc->sqrtsd(result.getXmm(), var.getMem());
		return result;
	}

	// Negation
	std::complex<double> negCC(std::complex<double>* args)
	{
		return -args[0];
	}

	template<>
	MpOperationNeg<double>::MpOperationNeg() :
		MpOperationFunc<double, double>(MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpFlagIsOperator, 1, VPTR(negRR), 3)
	{
	}
	template<>
	MpOperationNeg<std::complex<double>>::MpOperationNeg() :
		MpOperationFunc<std::complex<double>, std::complex<double>>(MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpFlagIsOperator, 1, VPTR(negCC), 3)
	{
	}

	template<typename T>
	JitVar MpOperationNeg<T>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->pxor(var.getXmm(), jc->getConstantU64(uint64_t(0x8000000000000000), uint64_t(0x8000000000000000)).getMem());
		return var;
	}

	template<typename T>
	uint32_t MpOperationNeg<T>::optimize(AstOptimizer * opt, AstNode * node) const
	{
		// as the reference to node might be invalidated by the call to MpOperationFunc::optimize,
		// we get the parent and the index of node within parent->_children.
		auto parent = node->getParent();
		size_t i;
		for (i = 0; i < parent->getLength(); i++)
		{
			if (parent->getAt(i) == node)
				break;
		}

		auto ret = MpOperationFunc<T, T>::optimize(opt, node);
		if (ret != ErrorCode::kErrorOk)
			return ret;

		// correct the reference to node.
		node = parent->getAt(i);

		// -(-(x)) = x
		if (node->getNodeType() == AstNodeType::kAstNodeUnaryOp &&
			node->getAt(0)->getNodeType() == AstNodeType::kAstNodeUnaryOp
			&& static_cast<AstUnaryOp*>(node)->_mpOp == static_cast<AstUnaryOp*>(node->getAt(0))->_mpOp)
		{
			AstNode* childOfChild = static_cast<AstUnaryOp*>(node->getAt(0))->unlinkChild();
			parent->replaceNode(node, childOfChild);
			opt->getAst()->deleteNode(node);
		}
		return ErrorCode::kErrorOk;
	}

	// Not
	std::complex<double> notCC(std::complex<double>* args)
	{
		return std::complex<double>(args[0] == std::complex<double>(0, 0) ? 1.0 : 0.0, 0.0);
	}

	template<>
	MpOperationNot<double>::MpOperationNot() :
		MpOperationFunc<double, double>(MpOperationFlags::OpFlagIsOperator, 1, VPTR(notRR), 3)
	{
	}
	template<>
	MpOperationNot<std::complex<double>>::MpOperationNot() :
		MpOperationFunc<std::complex<double>, std::complex<double>>(MpOperationFlags::OpFlagIsOperator, 1, VPTR(notCC), 3)
	{
	}

	template<>
	JitVar MpOperationNot<double>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVar(var);
		jc->cc->cmpsd(var.getXmm(), jc->getConstantD64AsPD(0.0).getMem(), int(asmjit::x86::kCmpEQ));
		jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return var;
	}

	template<>
	JitVar MpOperationNot<std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->cmppd(var.getXmm(), jc->getConstantD64(std::complex<double>(0.0, 0.0)).getMem(), int(asmjit::x86::kCmpEQ));
		jc->cc->andpd(var.getXmm(), jc->getConstantD64(std::complex<double>(1.0, 0.0)).getMem());
		return var;
	}

	// Conjugate
	std::complex<double> conjugCC(std::complex<double>* args) { return std::complex<double>(args->real(), -args->imag()); }

	MpOperationConjug::MpOperationConjug() : MpOperationFunc<std::complex<double>, std::complex<double>>(MpOperationFlags::OpFlagNone, 1, VPTR(conjugCC))
	{
	}

	JitVar MpOperationConjug::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar tmp = jc->onNode(node->getAt(0));
		JitVar result = jc->registerVarComplex(tmp, !node->getAt(0)->returnsComplex());
		jc->cc->pxor(result.getXmm(), jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000)).getMem());
		return result;
	}

	uint32_t MpOperationConjug::optimize(AstOptimizer * opt, AstNode * node) const
	{
		// as the reference to node might be invalidated by the call to MpOperationFunc::optimize,
		// we get the parent and the index of node within parent->_children.
		auto parent = node->getParent();
		size_t i;
		for (i = 0; i < parent->getLength(); i++)
		{
			if (parent->getAt(i) == node)
				break;
		}

		auto ret = MpOperationFunc<std::complex<double>, std::complex<double>>::optimize(opt, node);
		if (ret != ErrorCode::kErrorOk)
			return ret;

		// correct the reference to node.
		node = parent->getAt(i);

		// conj(conj(x)) = x
		if (node->getNodeType() == AstNodeType::kAstNodeUnaryOp &&
			node->getAt(0)->getNodeType() == AstNodeType::kAstNodeUnaryOp
			&& static_cast<AstUnaryOp*>(node)->_mpOp == static_cast<AstUnaryOp*>(node->getAt(0))->_mpOp)
		{
			AstNode* childOfChild = static_cast<AstUnaryOp*>(node->getAt(0))->unlinkChild();
			parent->replaceNode(node, childOfChild);
			opt->getAst()->deleteNode(node);
		}
		return ErrorCode::kErrorOk;
	}

	// Reciprocal
	std::complex<double> recipCC(std::complex<double>* args)
	{
		return 1.0 / args[0];
	}

	template<>
	MpOperationRecip<double>::MpOperationRecip() : MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(recipRR))
	{
	}
	template<>
	MpOperationRecip<std::complex<double>>::MpOperationRecip() : MpOperationFunc<std::complex<double>, std::complex<double>>(MpOperationFlags::OpFlagNone, 1, VPTR(recipCC))
	{
	}

	template<>
	JitVar MpOperationRecip<double>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		JitVar result;
		result = JitVar(jc->cc->newXmmSd(), false);
		jc->cc->movsd(result.getXmm(), jc->getConstantD64(1.0).getMem());
		if (var.isMem())
			jc->cc->divsd(result.getXmm(), var.getMem());
		else
			jc->cc->divsd(result.getXmm(), var.getXmm());
		return result;
	}
	template<>
	JitVar MpOperationRecip<std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		JitVar result;
		// as of http://www.chemistrylearning.com/reciprocal-of-a-complex-number/
		var = jc->writableVarComplex(var);
		result = JitVar(jc->cc->newXmmPd(), false);
		jc->cc->movapd(result.getXmm(), var.getXmm());
		jc->cc->mulpd(var.getXmm(), var.getXmm());
		jc->cc->haddpd(var.getXmm(), var.getXmm());
		jc->cc->pxor(result.getXmm(), jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000)).getMem());
		jc->cc->divpd(result.getXmm(), var.getXmm());
		return result;
	}

	// sign bit
	MpOperationSignBit::MpOperationSignBit() : MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(signbitRR))
	{
	}

	JitVar MpOperationSignBit::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar result(jc->cc->newXmmSd(), false);
		jc->cc->pshufd(result.getXmm(), jc->registerVar(var).getXmm(), asmjit::x86::shufImm(3, 2, 1, 1));
		jc->cc->psrad(result.getXmm(), 31);
		jc->cc->andpd(result.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return result;
	}

	// Copy sign
	MpOperationCopySign::MpOperationCopySign() : MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 2, VPTR(copysignRR))
	{
	}

	JitVar MpOperationCopySign::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar vl = jc->writableVar(jc->onNode(node->getAt(0)));
		JitVar vr = jc->writableVar(jc->onNode(node->getAt(1)));
		jc->cc->andpd(vl.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x7FFFFFFFFFFFFFFF)).getMem());
		jc->cc->andpd(vr.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->orpd(vl.getXmm(), vr.getXmm());

		return vl;
	}


	// Average
	std::complex<double> avgCC(std::complex<double>* args)
	{
		return (args[0] + args[1]) * 0.5;
	}

	template<>
	MpOperationAvg<double>::MpOperationAvg() : MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 2, VPTR(avgRR))
	{
	}
	template<>
	MpOperationAvg<std::complex<double>>::MpOperationAvg() : MpOperationFunc<std::complex<double>, std::complex<double>>(MpOperationFlags::OpFlagNone, 2, VPTR(avgCC))
	{
	}

	template<>
	JitVar MpOperationAvg<double>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar vl = jc->onNode(node->getAt(0));;
		JitVar vr = jc->onNode(node->getAt(1));

		vl = jc->writableVar(jc->onNode(node->getAt(0)));
		vr = jc->onNode(node->getAt(1));
		if (vr.isMem())
		{
			jc->cc->addsd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->addsd(vl.getXmm(), vr.getXmm());
		}
		jc->cc->mulsd(vl.getXmm(), jc->getConstantD64(0.5).getXmm());
		return vl;
	}
	template<>
	JitVar MpOperationAvg<std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar vl = jc->onNode(node->getAt(0));;
		JitVar vr = jc->onNode(node->getAt(1));

		if (!node->getAt(0)->returnsComplex())
		{
			vl = jc->registerVarAsComplex(vl);
		}
		else
		{
			vl = jc->writableVarComplex(vl);
		}

		if (!node->getAt(1)->returnsComplex())
		{
			vr = jc->registerVarAsComplex(vr);
		}

		if (vr.isMem())
		{
			jc->cc->addpd(vl.getXmm(), vr.getMem());
		}
		else
		{
			jc->cc->addpd(vl.getXmm(), vr.getXmm());
		}
		jc->cc->mulpd(vl.getXmm(), jc->getConstantD64(std::complex<double>(0.5, 0.5)).getXmm());
		return vl;
	}


	// Absolute
	MpOperationAbs::MpOperationAbs() : MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(absRR))
	{
	}

	JitVar MpOperationAbs::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar result;
		var = jc->writableVar(var);
		result = JitVar(jc->cc->newXmmSd(), false);
		jc->cc->xorpd(result.getXmm(), result.getXmm());
		jc->cc->subsd(result.getXmm(), var.getXmm());
		jc->cc->maxsd(result.getXmm(), var.getXmm());
		return result;
	}


	// round
	MpOperationRound::MpOperationRound() : MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(roundRR))
	{
	}

	JitVar MpOperationRound::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), false);

		if (jc->enableSSE4_1)
		{
			JitVar tmp(jc->cc->newXmmSd(), false);
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

			JitVar t1(jc->cc->newXmmSd(), false);
			JitVar t2(jc->cc->newXmmSd(), false);
			JitVar t3(jc->cc->newXmmSd(), false);

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
	MpOperationRoundEven::MpOperationRoundEven() :MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(roundevenRR))
	{
	}

	JitVar MpOperationRoundEven::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), false);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(result.getXmm(), var.getXmm(), asmjit::x86::kRoundNearest | asmjit::x86::kRoundInexact);
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), false);
			JitVar t2(jc->cc->newXmmSd(), false);

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
	MpOperationTrunc::MpOperationTrunc() :MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(truncRR))
	{
	}

	JitVar MpOperationTrunc::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), false);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(result.getXmm(), var.getXmm(), asmjit::x86::kRoundTrunc | asmjit::x86::kRoundInexact);
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), false);
			JitVar t2(jc->cc->newXmmSd(), false);
			JitVar t3(jc->cc->newXmmSd(), false);

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
	MpOperationFrac::MpOperationFrac() :MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(fracRR))
	{
	}

	JitVar MpOperationFrac::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar tmp(jc->cc->newXmmSd(), false);

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

			JitVar t1(jc->cc->newXmmSd(), false);
			JitVar t2(jc->cc->newXmmSd(), false);
			JitVar t3(jc->cc->newXmmSd(), false);

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
	MpOperationFloor::MpOperationFloor() :MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(floorRR))
	{
	}

	JitVar MpOperationFloor::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), false);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(result.getXmm(), var.getXmm(), asmjit::x86::kRoundDown | asmjit::x86::kRoundInexact);
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), false);
			JitVar t2(jc->cc->newXmmSd(), false);
			JitVar t3(jc->cc->newXmmSd(), false);

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
	MpOperationcCeil::MpOperationcCeil() : MpOperationFunc<double, double>(MpOperationFlags::OpFlagNone, 1, VPTR(ceilRR))
	{
	}

	JitVar MpOperationcCeil::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar var(jc->writableVar(jc->onNode(node->getAt(0))));
		JitVar result(jc->cc->newXmmSd(), false);

		if (jc->enableSSE4_1)
		{
			jc->cc->roundsd(result.getXmm(), var.getXmm(), asmjit::x86::kRoundUp | asmjit::x86::kRoundInexact);
		}
		else
		{
			const double maxn = 4503599627370496.0;
			const double magic0 = 6755399441055744.0;

			JitVar t1(jc->cc->newXmmSd(), false);
			JitVar t2(jc->cc->newXmmSd(), false);
			JitVar t3(jc->cc->newXmmSd(), false);

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

	template<>
	JitVar MpOperationBinary<double>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar vl, vr;
		AstNode* left = node->getAt(0);
		AstNode* right = node->getAt(1);

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
		return generateAsm(jc, vl, vr);
	}

	template<>
	JitVar MpOperationBinary<std::complex<double>>::compile(JitCompiler * jc, AstNode * node) const
	{
		JitVar vl, vr;
		AstNode* left = node->getAt(0);
		AstNode* right = node->getAt(1);

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
		return generateAsm(jc, vl, vr);
	}

	template<typename T>
	uint32_t MpOperationBinary<T>::optimize(AstOptimizer * opt, AstNode * node) const
	{
		AstNode* left = node->getAt(0);
		AstNode* right = node->getAt(1);

		bool lIsImm = left->isImm();
		bool rIsImm = right->isImm();

		if (lIsImm && rIsImm && !hasFlag(MpOperationFlags::OpFlagHasState))
		{
			// optimize a calculation with two immediates.
			AstImm* lNode = static_cast<AstImm*>(left);
			AstImm* rNode = static_cast<AstImm*>(right);

			lNode->setValue(calculate(lNode->getValue<T>(), rNode->getValue<T>()));

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
			if ((hasFlag(MpOperationFlags::OpFlagNopIfLZero) && lNode->getValue<std::complex<double>>() == std::complex<double>(0.0, 0.0)) ||
				(hasFlag(MpOperationFlags::OpFlagNopIfLOne) && lNode->getValue<std::complex<double>>() == std::complex<double>(1.0, 0.0)))
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

			if ((hasFlag(MpOperationFlags::OpFlagNopIfRZero) && rNode->getValue<std::complex<double>>() == std::complex<double>(0.0, 0.0)) ||
				(hasFlag(MpOperationFlags::OpFlagNopIfROne) && rNode->getValue<std::complex<double>>() == std::complex<double>(1.0, 0.0)))
			{
				node->_children[0]->_parent = nullptr;
				node->_children[0] = nullptr;
				node->getParent()->replaceNode(node, left);
				opt->getAst()->deleteNode(node);
			}

		}
		return ErrorCode::kErrorOk;
	}

	template<typename T>
	JitVar MpOperationBinary<T>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
	{
		throw std::runtime_error("No Override available!");
	}

	template<typename T>
	T MpOperationBinary<T>::calculate(T vl, T vr) const
	{
		return std::numeric_limits<T>::quiet_NaN();
	}

	template<>
	std::complex<double> MpOperationBinary<std::complex<double>>::calculate(std::complex<double> vl, std::complex<double> vr)  const
	{
		return std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
	}

	template<>
	MpOperationAdd<double>::MpOperationAdd() :
		MpOperationBinary<double>(Signature(2, MpOperationFlags::OpFlagNopIfZero | MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagIsOperator), 6)
	{
	}
	template<>
	MpOperationAdd<std::complex<double>>::MpOperationAdd() :
		MpOperationBinary<std::complex<double>>(Signature(2, Signature::type::complex, MpOperationFlags::OpFlagNopIfZero | MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagIsOperator), 6)
	{
	}

	// Addition
	template<>
	JitVar MpOperationAdd<double>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	template<>
	JitVar MpOperationAdd<std::complex<double>>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	template<typename T>
	T MpOperationAdd<T>::calculate(T vl, T vr) const
	{
		return vl + vr;
	}


	// Subtraction
	template<>
	MpOperationSub<double>::MpOperationSub() :
		MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpFlagNopIfRZero | MpOperationFlags::OpFlagIsOperator), 6)
	{
	}
	template<>
	MpOperationSub<std::complex<double>>::MpOperationSub() :
		MpOperationBinary<std::complex<double>>(Signature(2, Signature::type::complex, MpOperationFlags::OpFlagNopIfRZero | MpOperationFlags::OpFlagIsOperator), 6)
	{
	}

	template<>
	JitVar MpOperationSub<double>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	template<>
	JitVar MpOperationSub<std::complex<double>>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	template<typename T>
	T MpOperationSub<T>::calculate(T vl, T vr) const
	{
		return vl - vr;
	}


	// Multiplication
	template<>
	MpOperationMul<double>::MpOperationMul() :
		MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpFlagNopIfZero | MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagIsOperator), 5)
	{
	}
	template<>
	MpOperationMul<std::complex<double>>::MpOperationMul() :
		MpOperationBinary<std::complex<double>>(Signature(2, Signature::type::complex, MpOperationFlags::OpFlagNopIfZero | MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagIsOperator), 5)
	{
	}

	template<>
	JitVar MpOperationMul<double>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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
	template<>
	JitVar MpOperationMul<std::complex<double>>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
	{
		if (vr.isMem())
		{
			vr = jc->writableVarComplex(vr);
		}
		if (vl == vr)
		{
			vr = jc->copyVarComplex(vl, false);
		}
		JitVar ret(jc->cc->newXmmPd(), false);
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

	template<typename T>
	T MpOperationMul<T>::calculate(T vl, T vr) const
	{
		return vl * vr;
	}

	// Division
	template<>
	MpOperationDiv<double>::MpOperationDiv() :
		MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpFlagNopIfLOne | MpOperationFlags::OpFlagIsOperator), 5)
	{
	}
	template<>
	MpOperationDiv<std::complex<double>>::MpOperationDiv() :
		MpOperationBinary<std::complex<double>>(Signature(2, Signature::type::complex, MpOperationFlags::OpFlagNopIfLOne | MpOperationFlags::OpFlagIsOperator), 5)
	{
	}
	template<>
	JitVar MpOperationDiv<double>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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
	template<>
	JitVar MpOperationDiv<std::complex<double>>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
	{
		if (vr.isMem())
		{
			vr = jc->writableVarComplex(vr);
		}
		if (vl == vr)
		{
			vr = jc->copyVarComplex(vl, false);
		}
		JitVar ret(jc->cc->newXmmPd(), false);
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

	template<typename T>
	T MpOperationDiv<T>::calculate(T vl, T vr) const
	{
		return vl / vr;
	}

	// Minimum
	JitVar MpOperationMin::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	double MpOperationMin::calculate(double vl, double vr)  const
	{
		return std::min(vl, vr);
	}

	// Maximum
	JitVar MpOperationMax::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	double MpOperationMax::calculate(double vl, double vr)  const
	{
		return std::max(vl, vr);
	}

	// Equality
	template<>
	MpOperationEq<double>::MpOperationEq() :
		MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagIsOperator), 9)
	{
	}
	template<>
	MpOperationEq<std::complex<double>>::MpOperationEq() :
		MpOperationBinary<std::complex<double>>(Signature(2, Signature::type::complex, MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagIsOperator), 9)
	{
	}

	template<>
	JitVar MpOperationEq<double>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	template<>
	JitVar MpOperationEq<std::complex<double>>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	template<typename T>
	T MpOperationEq<T>::calculate(T vl, T vr)  const
	{
		return vl == vr ? 1.0 : 0.0;
	}

	// Inequality
	template<>
	MpOperationNe<double>::MpOperationNe() :
		MpOperationBinary<double>(Signature(2, Signature::type::real, MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagIsOperator), 9)
	{
	}
	template<>
	MpOperationNe<std::complex<double>>::MpOperationNe() :
		MpOperationBinary<std::complex<double>>(Signature(2, Signature::type::complex, MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagIsOperator), 9)
	{
	}

	template<>
	JitVar MpOperationNe<double>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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
	template<>
	JitVar MpOperationNe<std::complex<double>>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	template<typename T>
	T MpOperationNe<T>::calculate(T vl, T vr)  const
	{
		return vl != vr ? 1.0 : 0.0;
	}

	// Lesser than
	JitVar MpOperationLt::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	double MpOperationLt::calculate(double vl, double vr)  const
	{
		return vl < vr ? 1.0 : 0.0;
	}

	// Lesser equal
	JitVar MpOperationLe::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	double MpOperationLe::calculate(double vl, double vr)  const
	{
		return vl <= vr ? 1.0 : 0.0;
	}

	// Greater than
	JitVar MpOperationGt::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	double MpOperationGt::calculate(double vl, double vr)  const
	{
		return vl > vr ? 1.0 : 0.0;
	}

	// Greater equal
	JitVar MpOperationGe::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	double MpOperationGe::calculate(double vl, double vr)  const
	{
		return vl >= vr ? 1.0 : 0.0;
	}

	// Modulo
	JitVar MpOperationModulo::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
	{
		JitVar result(jc->cc->newXmmSd(), false);
		JitVar tmp(jc->cc->newXmmSd(), false);

		vl = jc->writableVar(vl);
		if (vl == vr)
		{
			vr = jc->copyVar(vl, false);
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

			JitVar t1(jc->cc->newXmmSd(), false);
			JitVar t2(jc->cc->newXmmSd(), false);
			JitVar t3(jc->cc->newXmmSd(), false);

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

	double MpOperationModulo::calculate(double vl, double vr) const
	{
		return fmod(vl, vr);
	}

	// Ternary operation

	template<>
	MpOperationTernary<double>::MpOperationTernary() : MpOperation(Signature(3, Signature::type::real, MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpFlagIsOperator), 15)
	{
	}
	template<>
	MpOperationTernary<std::complex<double>>::MpOperationTernary() : MpOperation(Signature(3, Signature::type::complex, MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpFlagIsOperator), 15)
	{
	}

	template<typename T>
	JitVar MpOperationTernary<T>::compile(JitCompiler* jc, AstNode * node) const
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
		bool lIsVarOrImm = left->getNodeType() == AstNodeType::kAstNodeVar || left->getNodeType() == AstNodeType::kAstNodeImm;

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
		bool rIsVarOrImm = right->getNodeType() == AstNodeType::kAstNodeVar || right->getNodeType() == AstNodeType::kAstNodeImm;

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
		if (node->hasNodeFlag(AstNodeFlags::kAstReturnsComplex))
			return jc->copyVarComplex(JitVar(regErg, false), false);
		else
			return jc->copyVar(JitVar(regErg, false), false);
	}

	template<typename T>
	uint32_t MpOperationTernary<T>::optimize(AstOptimizer *opt, AstNode *node) const
	{
		AstTernaryOp * ternaryNode = reinterpret_cast<AstTernaryOp*>(node);
		AstNode* branchCond = ternaryNode->getCondition();
		if (branchCond->isImm())
		{
			// optimize an immediate condition
			bool conditionIsTrue = static_cast<AstImm*>(branchCond)->getValue<std::complex<double>>() != std::complex<double>({ 0, 0 });

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

			opt->_ast->deleteNode(ternaryNode);
		}

		return ErrorCode::kErrorOk;

	}

	// Assignment
	template<>
	MpOperationAssignment<double>::MpOperationAssignment() : MpOperation(Signature(1, Signature::type::real, MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpIsAssgignment), 15)
	{
	}
	template<>
	MpOperationAssignment<std::complex<double>>::MpOperationAssignment() : MpOperation(Signature(1, Signature::type::complex, MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpIsAssgignment), 15)
	{
	}

	template<typename T>
	JitVar MpOperationAssignment<T>::compile(JitCompiler * jc, AstNode * node) const
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

	template<typename T>
	uint32_t MpOperationAssignment<T>::optimize(AstOptimizer * opt, AstNode * node) const
	{
		AstVarDecl * varDecl;
		if (node->getNodeType() == AstNodeType::kAstNodeVarDecl)
		{
			varDecl = static_cast<AstVarDecl*>(node);
		}
		else
		{
			return ErrorCode::kErrorInvalidState;
		}

		AstSymbol* sym = varDecl->getSymbol();


		if (varDecl->hasChild())
		{
			AstNode* child = varDecl->getChild();

			if (child->returnsComplex())
			{
				sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsComplex);
			}

			if (child->isImm())
			{
				sym->setValue(static_cast<AstImm*>(child)->getValue<T>());

				sym->setAssigned();
			}
		}

		return ErrorCode::kErrorOk;
	}

} // end namespace mathpresso
