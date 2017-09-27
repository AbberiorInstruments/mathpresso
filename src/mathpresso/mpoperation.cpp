// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

#include <mathpresso/mpoperation.h>
#include <mathpresso/mpast_p.h>
#include <mathpresso/mpcompiler_p.h>
#include <mathpresso/mpeval_p.h>
#include <mathpresso/mpoptimizer_p.h>
#include <asmjit/x86/x86operand.h>
#include <asmjit/x86/x86inst.h>

// Missing in std
namespace std
{
	cplx_t log2(const cplx_t &x) { return std::log(x) / std::log(2.); }
}

#define _ADD_FUNRC(name, fptr) {\
	static auto p = static_cast<double(*)(const cplx_t&)>(fptr); \
	ctx ->addObject(name, _OBJ(static_cast<double(*)(const cplx_t&)>(fptr))); }

#define _ADD_CPLX1(name, fptr) {\
	static auto p = static_cast<cplx_t(*)(const cplx_t&)>(fptr); \
	ctx ->addObject(name, _OBJ(static_cast<cplx_t(*)(const cplx_t&)>(fptr))); }

#define _ADD_CPLX2(name, fptr) {\
	static auto p = static_cast<cplx_t(*)(const cplx_t&, const cplx_t&)>(fptr); \
	ctx ->addObject(name, _OBJ(static_cast<cplx_t(*)(const cplx_t&, const cplx_t&)>(fptr))); }

// First line is to ensure the instantiation of the function
#define _ADD_FUNC1(name, fptr) _ADD_CPLX1(name, fptr) \
	ctx ->addObject(name, _OBJ(static_cast<double(*)(double)>(fptr))); 

#define _ADD_FUNC2(name, fptr) _ADD_CPLX2(name, fptr) \
	ctx ->addObject(name, _OBJ(static_cast<double(*)(double, double)>(fptr)));

namespace mathpresso
{
	std::complex<double> sqrtc(double x)
	{
		return std::sqrt(std::complex<double>(x, 0));
	}

	Signature::Signature(size_t nargs, type cmnType) noexcept
	{
		init(cmnType, { nargs, { cmnType, "" } });
	}

	Signature::Signature(size_t nargs, type argType, type retType) noexcept
	{
		init(retType, { nargs,{ argType, "" } });
	}

	Signature::Signature(size_t nargs) noexcept
	{
		init(type::real, { nargs, { type::real, "" } });
	}

	bool Signature::areParams(type _type) const
	{
		bool ret = true;
		for (auto p : parameters_)
		{
			ret &= (p.type_ == _type);
		}
		return ret;
	}

	std::string Signature::to_string()
	{
		std::string out("");
		out += typeToString(return_type_);
		out += " (";
		for (size_t i = 0; i < parameters_.size(); i++)
		{
			out += typeToString(parameters_[i].type_);
			if (parameters_[i].name_ != "")
				out += " " + parameters_[i].name_;
			if (i != parameters_.size() - 1)
				out += ", ";
		}
		out += ")";
		return out;
	}

	std::string Signature::typeToString(type _type)
	{
		switch (_type)
		{
			case type::real: return "real";
			case type::complex: return "complex";
			default:
				throw std::runtime_error("unknown type.");
		}
	}

	void Signature::init(type retType, std::vector<param> params)
	{
		return_type_ = retType;
		parameters_ = params;
	}


	template<>
	MpOperationEval<double, double>::MpOperationEval(size_t numargs, uint32_t flags, uint32_t priority)  noexcept
		: MpOperation(Signature(numargs, Signature::type::real), flags, priority)
	{
	}

	template<>
	MpOperationEval<cplx_t, double>::MpOperationEval(size_t numargs, uint32_t flags, uint32_t priority) noexcept
		: MpOperation(Signature(numargs, Signature::type::real, Signature::type::complex), flags, priority)
	{
	}

	template<>
	MpOperationEval<double, cplx_t>::MpOperationEval(size_t numargs, uint32_t flags, uint32_t priority) noexcept
		: MpOperation(Signature(numargs, Signature::type::complex, Signature::type::real), flags, priority)
	{
	}

	template<>
	MpOperationEval<cplx_t, cplx_t>::MpOperationEval(size_t numargs, uint32_t flags, uint32_t priority) noexcept
		: MpOperation(Signature(numargs, Signature::type::complex), flags, priority)
	{
	}


	template<>
	JitVar MpOperationFunc<double, double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		asmjit::X86Xmm result = jc->cc->newXmmSd();
		asmjit::X86Xmm args[8];
		if (!fnPtr_)
		{
			// Should never happen, as the optimizer should have taken care of that. Remove later
			throw std::runtime_error("Implementation error!");
		}

		for (size_t i = 0; i < nargs(); i++)
		{
			args[i] = jc->registerVar(jc->onNode(node->getAt(i))).getXmm();
		}
		jc->inlineCall<double, double>(result, args, nargs(), fnPtr_);

		return JitVar(result, false);
	}
	template<>
	JitVar MpOperationFunc<double, cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		asmjit::X86Xmm result = jc->cc->newXmmSd();
		asmjit::X86Xmm args[8];
		if (!fnPtr_)
		{
			// Should never happen, as the optimizer should have taken care of that. Remove later
			throw std::runtime_error("Implementation error!");
		}

		for (size_t i = 0; i < nargs(); i++)
		{
			args[i] = jc->registerVarComplex(jc->onNode(node->getAt(i)), !node->getAt(i)->returnsComplex()).getXmm();
		}

		jc->inlineCall<double, cplx_t>(result, args, nargs(), fnPtr_);

		return JitVar(result, false);
	}
	template<>
	JitVar MpOperationFunc<cplx_t, double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		asmjit::X86Xmm result = jc->cc->newXmmPd();
		asmjit::X86Xmm args[8];
		if (!fnPtr_)
		{
			// Should never happen, as the optimizer should have taken care of that. Remove later
			throw std::runtime_error("Implementation error!");
		}

		for (size_t i = 0; i < nargs(); i++)
		{
			args[i] = jc->registerVar(jc->onNode(node->getAt(i))).getXmm();
		}
		jc->inlineCall<cplx_t, double>(result, args, nargs(), fnPtr_);

		return JitVar(result, false);
	}
	template<>
	JitVar MpOperationFunc<cplx_t, cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		asmjit::X86Xmm result = jc->cc->newXmmPd();
		asmjit::X86Xmm args[8];
		if (!fnPtr_)
		{
			// Should never happen, as the optimizer should have taken care of that. Remove later
			throw std::runtime_error("Implementation error!");
		}

		for (size_t i = 0; i < nargs(); i++)
		{
			args[i] = jc->registerVarComplex(jc->onNode(node->getAt(i)), !node->getAt(i)->returnsComplex()).getXmm();
		}

		jc->inlineCall<cplx_t, cplx_t>(result, args, nargs(), fnPtr_);

		return JitVar(result, false);
	}

	template<typename RET, typename PARAM>
	uint32_t MpOperationEval<RET, PARAM>::optimize(AstOptimizer * opt, std::shared_ptr<AstNode> node) const
	{
		bool b_all_imm = true;

		// Gather Information about the child-nodes.
		for (size_t i = 0; i < node->getLength(); i++)
		{
			b_all_imm &= node->getAt(i)->isImm();
		}

		// optimize all-immediate calls:
		if (b_all_imm && !hasFlag(MpOperation::FlagHasState))
		{
			std::shared_ptr<AstImm> ret = std::make_shared<AstImm>(0);

			std::vector<PARAM> args;
			for (size_t i = 0; i < nargs(); i++)
			{
				args.push_back((std::static_pointer_cast<AstImm>(node->getAt(i)))->getValue<PARAM>());
			}

			ret->setValue(evaluate(args.data()));

			node->getParent()->replaceNode(node, ret);
			node = ret;
		}
		return ErrorCode::kErrorOk;
	}

	template<>
	double MpOperationFunc<double, double>::evaluate(const double * args) const
	{
		if (!fnPtr_)
		{
			throw std::runtime_error("Function does not exist.");
		}
#ifndef MATHPRESSO_ORIGINAL_DOUBLE_FUNCTION_CALLS
		return ((mpFuncDtoD)fnPtr_)(args);
#else
		switch (nargs())
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
	cplx_t MpOperationFunc<cplx_t, double>::evaluate(const double * args) const
	{
		if (!fnPtr_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpDtoC)fnPtr_)(args);
	}
	template<>
	double MpOperationFunc<double, cplx_t>::evaluate(const cplx_t* args) const
	{
		if (!fnPtr_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpCtoD)fnPtr_)(args);
	}
	template<>
	cplx_t MpOperationFunc<cplx_t, cplx_t>::evaluate(const cplx_t* args) const
	{
		if (!fnPtr_)
		{
			throw std::runtime_error("Function does not exist.");
		}
		return ((mpFuncpCtoC)fnPtr_)(args);
	}

	template<typename T>
	class MpOperationIsFinite : public MpOperationEval<double, T>
	{
	public:
		MpOperationIsFinite() : MpOperationEval<double, T>(1)
		{
		}
		virtual double evaluate(const T * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	template<>
	double MpOperationIsFinite<cplx_t>::evaluate(const cplx_t * args) const
	{
		return std::isfinite(args[0].real()) ? 1.0 : 0.0 && std::isfinite(args[0].imag()) ? 1.0 : 0.0;
	}

	template<>
	double MpOperationIsFinite<double>::evaluate(const double * args) const
	{
		return std::isfinite(args[0]) ? 1.0 : 0.0;
	}

	template<>
	JitVar MpOperationIsFinite<double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVar(var);
		jc->cc->orpd(var.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->cmpsd(var.getXmm(), jc->getConstantD64(0.0).getMem(), asmjit::x86::kCmpLE);
		jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return var;
	}

	template<>
	JitVar MpOperationIsFinite<cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->orpd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0x8000000000000000), MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->cmppd(var.getXmm(), jc->getConstantD64(cplx_t(0.0, 0.0)).getMem(), asmjit::x86::kCmpLE);
		jc->cc->andpd(var.getXmm(), jc->getConstantD64(cplx_t(1.0, 1.0)).getMem());

		JitVar tmp(jc->cc->newXmmSd(), false);
		jc->cc->movhlps(tmp.getXmm(), var.getXmm());
		jc->cc->maxsd(tmp.getXmm(), var.getXmm());

		return tmp;
	}


	// MpOperationIsInFinite
	template<typename T>
	class MpOperationIsInfinite : public MpOperationEval<double, T>
	{
	public:
		MpOperationIsInfinite() noexcept : MpOperationEval<double, T>(1)
		{
		}
		virtual double evaluate(const T * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};


	template<>
	double MpOperationIsInfinite<double>::evaluate(const double * args) const
	{
		return std::isinf(args[0]) ? 1.0 : 0.0;
	}

	template<>
	double MpOperationIsInfinite<cplx_t>::evaluate(const cplx_t * args) const
	{
		return std::isinf(args[0].real()) ? 1.0 : 0.0 && std::isinf(args[0].imag()) ? 1.0 : 0.0;
	}

	template<>
	JitVar MpOperationIsInfinite<double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVar(var);
		jc->cc->orpd(var.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->cmpsd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0xFFF0000000000000)).getMem(), asmjit::x86::kCmpEQ);
		jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return var;
	}

	template<>
	JitVar MpOperationIsInfinite<cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->orpd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0x8000000000000000), MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->cmppd(var.getXmm(), jc->getConstantU64(MATHPRESSO_UINT64_C(0xFFF0000000000000), MATHPRESSO_UINT64_C(0xFFF0000000000000)).getMem(), asmjit::x86::kCmpEQ);
		jc->cc->andpd(var.getXmm(), jc->getConstantD64(cplx_t(1.0, 1.0)).getMem());

		JitVar tmp(jc->cc->newXmmSd(), false);
		jc->cc->movhlps(tmp.getXmm(), var.getXmm());
		jc->cc->maxsd(tmp.getXmm(), var.getXmm());

		return tmp;
	}


	// MpOperationIsNan	
	template<typename T>
	class MpOperationIsNan : public MpOperationEval<double, T>
	{
	public:
		MpOperationIsNan() noexcept : MpOperationEval<double, T>(1)
		{
		}
		virtual double evaluate(const T * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	template<>
	double MpOperationIsNan<double>::evaluate(const double * args) const
	{
		return std::isnan(args[0]) ? 1.0 : 0.0;
	}

	template<>
	double MpOperationIsNan<cplx_t>::evaluate(const cplx_t * args) const
	{
		return std::isnan(args[0].real()) ? 1.0 : 0.0 && std::isnan(args[0].imag()) ? 1.0 : 0.0;
	}

	template<>
	JitVar MpOperationIsNan<double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVar(var);
		jc->cc->cmpsd(var.getXmm(), var.getXmm(), asmjit::x86::kCmpEQ); // compare of NaN with NaN is false
		jc->cc->andnpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return var;
	}

	template<>
	JitVar MpOperationIsNan<cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->cmppd(var.getXmm(), var.getXmm(), asmjit::x86::kCmpEQ);
		jc->cc->andnpd(var.getXmm(), jc->getConstantD64(cplx_t(1.0, 1.0)).getMem());

		JitVar tmp(jc->cc->newXmmSd(), false);
		jc->cc->movhlps(tmp.getXmm(), var.getXmm());
		jc->cc->maxsd(tmp.getXmm(), var.getXmm());

		return tmp;
	}

	// MpOperationGetReal
	class MpOperationGetReal : public MpOperationEval<double, cplx_t>
	{
	public:
		MpOperationGetReal() noexcept : MpOperationEval<double, cplx_t>(1)
		{
		}
		virtual double evaluate(const cplx_t * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationGetReal::evaluate(const cplx_t * args) const
	{
		return args->real();
	}

	JitVar MpOperationGetReal::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	class MpOperationGetImag : public MpOperationEval<double, cplx_t>
	{
	public:
		MpOperationGetImag() noexcept : MpOperationEval<double, cplx_t>(1)
		{
		}
		virtual double evaluate(const cplx_t * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationGetImag::evaluate(const cplx_t * args) const
	{
		return args->imag();
	}

	JitVar MpOperationGetImag::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar varRet(jc->cc->newXmmSd(), false);;
		var = jc->registerVarComplex(var, !node->getAt(0)->hasNodeFlag(AstNodeFlags::kAstReturnsComplex));
		jc->cc->xorpd(varRet.getXmm(), varRet.getXmm());
		jc->cc->shufpd(var.getXmm(), varRet.getXmm(), asmjit::x86::shufImm(0, 1));
		return var;
	}

	// Square root
	class MpOperationSqrt : public MpOperationEval<double, double>
	{
	public:
		MpOperationSqrt() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationSqrt::evaluate(const double * args) const
	{
		return std::sqrt(args[0]);
	}

	JitVar MpOperationSqrt::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	template<typename T>
	class MpOperationNeg : public MpOperationEval<T, T>
	{
	public:
		MpOperationNeg() noexcept : MpOperationEval<T, T>(1)
		{
		}
		virtual T evaluate(const T * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	protected:
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};

	template<>
	double MpOperationNeg<double>::evaluate(const double * args) const
	{
		return -args[0];
	}

	template<>
	cplx_t MpOperationNeg<cplx_t>::evaluate(const cplx_t * args) const
	{
		return -args[0];
	}

	template<typename T>
	JitVar MpOperationNeg<T>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->pxor(var.getXmm(), jc->getConstantU64(uint64_t(0x8000000000000000), uint64_t(0x8000000000000000)).getMem());
		return var;
	}

	template<typename T>
	uint32_t MpOperationNeg<T>::optimize(AstOptimizer * opt, std::shared_ptr<AstNode> node) const
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

		auto ret = MpOperationEval<T, T>::optimize(opt, node);
		if (ret != ErrorCode::kErrorOk)
			return ret;

		// correct the reference to node.
		node = parent->getAt(i);

		// -(-(x)) = x
		if (node->getNodeType() == AstNodeType::kAstNodeUnaryOp &&
			node->getAt(0)->getNodeType() == AstNodeType::kAstNodeUnaryOp
			&& std::static_pointer_cast<AstUnaryOp>(node)->_mpOp == std::static_pointer_cast<AstUnaryOp>(node->getAt(0))->_mpOp)
		{
			std::shared_ptr<AstNode> childOfChild = std::static_pointer_cast<AstUnaryOp>(node->getAt(0))->unlinkChild();
			parent->replaceNode(node, childOfChild);
		}
		return ErrorCode::kErrorOk;
	}

	// Not
	template<typename T>
	class MpOperationNot : public MpOperationEval<T, T>
	{
	public:
		MpOperationNot() noexcept : MpOperationEval<T, T>(1)
		{
		}
		virtual T evaluate(const T * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	template<>
	double MpOperationNot<double>::evaluate(const double * args) const
	{
		return args[0] == 0.0 ? 1.0 : 0.0;
	}

	template<>
	cplx_t MpOperationNot<cplx_t>::evaluate(const cplx_t * args) const
	{
		return cplx_t(args[0] == cplx_t(0, 0) ? 1.0 : 0.0, 0.0);
	}

	template<>
	JitVar MpOperationNot<double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVar(var);
		jc->cc->cmpsd(var.getXmm(), jc->getConstantD64AsPD(0.0).getMem(), int(asmjit::x86::kCmpEQ));
		jc->cc->andpd(var.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return var;
	}

	template<>
	JitVar MpOperationNot<cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var = jc->onNode(node->getAt(0));
		var = jc->writableVarComplex(var);
		jc->cc->cmppd(var.getXmm(), jc->getConstantD64(cplx_t(0.0, 0.0)).getMem(), int(asmjit::x86::kCmpEQ));
		jc->cc->andpd(var.getXmm(), jc->getConstantD64(cplx_t(1.0, 0.0)).getMem());
		return var;
	}

	// Conjugate
	class MpOperationConjug : public MpOperationEval<cplx_t, cplx_t>
	{
	public:
		MpOperationConjug() noexcept : MpOperationEval<cplx_t, cplx_t>(1)
		{
		}
		virtual cplx_t evaluate(const cplx_t * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	protected:
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};

	cplx_t MpOperationConjug::evaluate(const cplx_t * args) const
	{
		return std::conj(args[0]);
	}

	JitVar MpOperationConjug::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar tmp = jc->onNode(node->getAt(0));
		JitVar result = jc->registerVarComplex(tmp, !node->getAt(0)->returnsComplex());
		jc->cc->pxor(result.getXmm(), jc->getConstantU64(uint64_t(0), uint64_t(0x8000000000000000)).getMem());
		return result;
	}

	uint32_t MpOperationConjug::optimize(AstOptimizer * opt, std::shared_ptr<AstNode> node) const
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

		auto ret = MpOperationEval<cplx_t, cplx_t>::optimize(opt, node);
		if (ret != ErrorCode::kErrorOk)
			return ret;

		// correct the reference to node.
		node = parent->getAt(i);

		// conj(conj(x)) = x
		if (node->getNodeType() == AstNodeType::kAstNodeUnaryOp &&
			node->getAt(0)->getNodeType() == AstNodeType::kAstNodeUnaryOp
			&& std::static_pointer_cast<AstUnaryOp>(node)->_mpOp == std::static_pointer_cast<AstUnaryOp>(node->getAt(0))->_mpOp)
		{
			std::shared_ptr<AstNode> childOfChild = std::static_pointer_cast<AstUnaryOp>(node->getAt(0))->unlinkChild();
			parent->replaceNode(node, childOfChild);
		}
		return ErrorCode::kErrorOk;
	}

	// Reciprocate
	template<typename T>
	class MpOperationRecip : public MpOperationEval<T, T>
	{
	public:
		MpOperationRecip() noexcept : MpOperationEval<T, T>(1)
		{
		}
		virtual T evaluate(const T * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	template<>
	double MpOperationRecip<double>::evaluate(const double * args) const
	{
		return 1.0 / args[0];
	}

	template<>
	cplx_t MpOperationRecip<cplx_t>::evaluate(const cplx_t * args) const
	{
		return 1.0 / args[0];
	}

	template<>
	JitVar MpOperationRecip<double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	JitVar MpOperationRecip<cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	class MpOperationSignBit : public MpOperationEval<double, double>
	{
	public:
		MpOperationSignBit() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double* args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationSignBit::evaluate(const double * args) const
	{
		return std::isinf(args[0]) ? 1.0 : 0.0;
	}

	JitVar MpOperationSignBit::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar var(jc->onNode(node->getAt(0)));
		JitVar result(jc->cc->newXmmSd(), false);
		jc->cc->pshufd(result.getXmm(), jc->registerVar(var).getXmm(), asmjit::x86::shufImm(3, 2, 1, 1));
		jc->cc->psrad(result.getXmm(), 31);
		jc->cc->andpd(result.getXmm(), jc->getConstantD64AsPD(1.0).getMem());
		return result;
	}

	// Copy sign
	class MpOperationCopySign : public MpOperationEval<double, double>
	{
	public:
		MpOperationCopySign() noexcept : MpOperationEval<double, double>(2)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationCopySign::evaluate(const double * args) const
	{
		return std::copysign(args[0], args[1]);
	}

	JitVar MpOperationCopySign::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar vl = jc->writableVar(jc->onNode(node->getAt(0)));
		JitVar vr = jc->writableVar(jc->onNode(node->getAt(1)));
		jc->cc->andpd(vl.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x7FFFFFFFFFFFFFFF)).getMem());
		jc->cc->andpd(vr.getXmm(), jc->getConstantU64AsPD(MATHPRESSO_UINT64_C(0x8000000000000000)).getMem());
		jc->cc->orpd(vl.getXmm(), vr.getXmm());

		return vl;
	}


	// Average
	template<typename T>
	class MpOperationAvg : public MpOperationEval<T, T>
	{
	public:
		MpOperationAvg() noexcept : MpOperationEval<T, T>(2)
		{
		}
		virtual T evaluate(const T * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	template<>
	double MpOperationAvg<double>::evaluate(const double * args) const
	{
		return (args[0] + args[1]) * 0.5;
	}

	template<>
	cplx_t MpOperationAvg<cplx_t>::evaluate(const cplx_t * args) const
	{
		return (args[0] + args[1]) * 0.5;
	}

	template<>
	JitVar MpOperationAvg<double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	JitVar MpOperationAvg<cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
		jc->cc->mulpd(vl.getXmm(), jc->getConstantD64(cplx_t(0.5, 0.5)).getXmm());
		return vl;
	}

	// Absolute
	class MpOperationAbs : public MpOperationEval<double, double>
	{
	public:
		MpOperationAbs() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationAbs::evaluate(const double * args) const
	{
		return std::abs(args[0]);
	}

	JitVar MpOperationAbs::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	class MpOperationRound : public MpOperationEval<double, double>
	{
	public:
		MpOperationRound() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationRound::evaluate(const double * args) const
	{
		return std::floor(args[0] + .5);
	}

	JitVar MpOperationRound::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	class MpOperationRoundEven : public MpOperationEval<double, double>
	{
	public:
		MpOperationRoundEven() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationRoundEven::evaluate(const double * args) const
	{
		return std::rint(args[0]);
	}

	JitVar MpOperationRoundEven::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	class MpOperationTrunc : public MpOperationEval<double, double>
	{
	public:
		MpOperationTrunc() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationTrunc::evaluate(const double * args) const
	{
		return std::trunc(args[0]);
	}

	JitVar MpOperationTrunc::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	class MpOperationFrac : public MpOperationEval<double, double>
	{
	public:
		MpOperationFrac() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationFrac::evaluate(const double * args) const
	{
		return args[0] - std::floor(args[0]);
	}

	JitVar MpOperationFrac::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	class MpOperationFloor : public MpOperationEval<double, double>
	{
	public:
		MpOperationFloor() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationFloor::evaluate(const double * args) const
	{
		return std::floor(args[0]);
	}

	JitVar MpOperationFloor::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	class MpOperationcCeil : public MpOperationEval<double, double>
	{
	public:
		MpOperationcCeil() noexcept : MpOperationEval<double, double>(1)
		{
		}
		virtual double evaluate(const double * args) const override;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
	};

	double MpOperationcCeil::evaluate(const double * args) const
	{
		return std::ceil(args[0]);
	}

	JitVar MpOperationcCeil::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
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
	JitVar MpOperationBinary<double>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar vl, vr;
		std::shared_ptr<AstNode> left = node->getAt(0);
		std::shared_ptr<AstNode> right = node->getAt(1);

		// check whether the vars are the same, to reduce memory-operations
		if (left->isVar() && right->isVar() &&
			std::static_pointer_cast<AstVar>(left)->getSymbol() == std::static_pointer_cast<AstVar>(right)->getSymbol())
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
	JitVar MpOperationBinary<cplx_t>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar vl, vr;
		std::shared_ptr<AstNode> left = node->getAt(0);
		std::shared_ptr<AstNode> right = node->getAt(1);

		// check whether the vars are the same, to reduce memory-operations
		if (left->isVar() && right->isVar() &&
			std::static_pointer_cast<AstVar>(left)->getSymbol() == std::static_pointer_cast<AstVar>(right)->getSymbol())
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
	uint32_t MpOperationBinary<T>::optimize(AstOptimizer * opt, std::shared_ptr<AstNode> node) const
	{
		std::shared_ptr<AstNode> left = node->getAt(0);
		std::shared_ptr<AstNode> right = node->getAt(1);

		bool lIsImm = left->isImm();
		bool rIsImm = right->isImm();

		if (lIsImm && rIsImm && !MpOperation::hasFlag(MpOperation::FlagHasState))
		{
			// optimize a calculation with two immediates.
			std::shared_ptr<AstImm> lNode = std::static_pointer_cast<AstImm>(left);
			std::shared_ptr<AstImm> rNode = std::static_pointer_cast<AstImm>(right);

			lNode->setValue(evaluate(lNode->getValue<T>(), rNode->getValue<T>()));

			// setValue sets the correct flags automatically.
			node->_children[0]->_parent.reset();
			node->_children[0] = nullptr;
			node->getParent()->replaceNode(node, lNode);
		}
		else if (lIsImm && !MpOperation::hasFlag(MpOperation::FlagHasState))
		{
			std::shared_ptr<AstImm> lNode = std::static_pointer_cast<AstImm>(left);
			// if the node is real, the imaginary part is set to zero by default.
			if ((MpOperation::hasFlag(NopIfLZero) && lNode->getValue<cplx_t>() == cplx_t(0.0, 0.0)) ||
				(MpOperation::hasFlag(NopIfLOne) && lNode->getValue<cplx_t>() == cplx_t(1.0, 0.0)))
			{
				node->_children[1]->_parent.reset();
				node->_children[1] = nullptr;
				node->getParent()->replaceNode(node, right);
			}
		}
		else if (rIsImm && !MpOperation::hasFlag(MpOperation::FlagHasState))
		{
			std::shared_ptr<AstImm> rNode = std::static_pointer_cast<AstImm>(right);

			if ((MpOperation::hasFlag(NopIfRZero) && rNode->getValue<cplx_t>() == cplx_t(0.0, 0.0)) ||
				(MpOperation::hasFlag(NopIfROne) && rNode->getValue<cplx_t>() == cplx_t(1.0, 0.0)))
			{
				node->_children[0]->_parent.reset();
				node->_children[0] = nullptr;
				node->getParent()->replaceNode(node, left);
			}

		}
		return ErrorCode::kErrorOk;
	}

	// Addition
	template<typename T>
	class MpOperationAdd : public MpOperationBinary<T>
	{
	public:
		MpOperationAdd() : MpOperationBinary<T>(MpOperationBinary<T>::NopIfZero | MpOperationBinary<T>::IsCommutativ, 6)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vL, T vR) const override
		{
			return vL + vR;
		}
	};

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
	JitVar MpOperationAdd<cplx_t>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	// Subtraction
	template<typename T>
	class MpOperationSub : public MpOperationBinary<T>
	{
	public:
		MpOperationSub() : MpOperationBinary<T>(MpOperationBinary<T>::NopIfZero, 6)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vL, T vR) const override
		{
			return vL - vR;
		}
	};

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
	JitVar MpOperationSub<cplx_t>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	// Multiplication
	template<typename T>
	class MpOperationMul : public MpOperationBinary<T>
	{
	public:
		MpOperationMul() : MpOperationBinary<T>(MpOperationBinary<T>::NopIfZero | MpOperationBinary<T>::IsCommutativ, 5)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vl, T vr) const override
		{
			return vl*vr;
		}
	};

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
	JitVar MpOperationMul<cplx_t>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	// Division
	template<typename T>
	class MpOperationDiv : public MpOperationBinary<T>
	{
	public:
		MpOperationDiv() : MpOperationBinary<T>(MpOperationBinary<T>::NopIfLOne, 5)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vl, T vr) const override
		{
			return vl / vr;
		}
	};

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
	JitVar MpOperationDiv<cplx_t>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	// Minimum
	class MpOperationMin : public MpOperationBinary<double>
	{
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return std::min(vl, vr);
		}
	};

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

	// Maximum
	class MpOperationMax : public  MpOperationBinary<double>
	{
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return std::max(vl, vr);
		}
	};

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

	// Equality
	template<typename T>
	class MpOperationEq : public MpOperationBinary<T>
	{
	public:
		MpOperationEq() : MpOperationBinary<T>(MpOperationBinary<T>::IsCommutativ, 9)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vl, T vr) const override
		{
			return vl == vr ? 1.0 : 0.0;
		}
	};

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
	JitVar MpOperationEq<cplx_t>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	// Inequality
	template<typename T>
	class MpOperationNe : public MpOperationBinary<T>
	{
	public:
		MpOperationNe() : MpOperationBinary<T>(MpOperationBinary<T>::IsCommutativ, 9)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual T evaluate(T vl, T vr) const override
		{
			return vl != vr ? 1.0 : 0.0;
		}
	};

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
	JitVar MpOperationNe<cplx_t>::generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const
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

	// Lesser than
	class MpOperationLt : public  MpOperationBinary<double>
	{
	public:
		MpOperationLt() : MpOperationBinary(MpOperation::None, 8)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return vl < vr ? 1. : 0.;
		}
	};

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

	// Lesser Equal
	class MpOperationLe : public  MpOperationBinary<double>
	{
	public:
		MpOperationLe() : MpOperationBinary<double>(MpOperation::None, 8)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return vl <= vr ? 1. : 0.;
		}
	};

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

	// Greater than
	class MpOperationGt : public  MpOperationBinary<double>
	{
	public:
		MpOperationGt() : MpOperationBinary<double>(MpOperation::None, 8)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return vl > vr ? 1. : 0.;
		}
	};

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


	// Greater equal
	class MpOperationGe : public  MpOperationBinary<double>
	{
	public:
		MpOperationGe() : MpOperationBinary<double>(MpOperation::None, 8)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return vl >= vr ? 1. : 0.;
		}
	};

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

	// Modulo
	class MpOperationModulo : public  MpOperationBinary<double>
	{
	public:
		MpOperationModulo() : MpOperationBinary<double>(MpOperation::None, 5)
		{
		}
	protected:
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double evaluate(double vl, double vr) const override
		{
			return std::fmod(vl, vr);
		}
	};
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

	// Ternary Operation
	template<typename T>
	class MpOperationTernary : public MpOperation
	{
	public:
		MpOperationTernary() noexcept;

		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	protected:
		bool isColon_;
	};

	template<>
	MpOperationTernary<double>::MpOperationTernary() noexcept
		: MpOperation(Signature(3, Signature::type::real), MpOperation::RighttoLeft, 15)
	{
	}
	template<>
	MpOperationTernary<cplx_t>::MpOperationTernary() noexcept
		: MpOperation(Signature(3, Signature::type::complex), MpOperation::RighttoLeft, 15)
	{
	}

	template<typename T>
	JitVar MpOperationTernary<T>::compile(JitCompiler* jc, std::shared_ptr<AstNode> node) const
	{
		asmjit::Label lblElse = jc->cc->newLabel();
		asmjit::Label lblEnd = jc->cc->newLabel();
		JitVar erg;
		std::shared_ptr<AstNode> left = std::static_pointer_cast<AstTernaryOp>(node)->getLeft();
		std::shared_ptr<AstNode> right = std::static_pointer_cast<AstTernaryOp>(node)->getRight();
		std::shared_ptr<AstNode> condition = std::static_pointer_cast<AstTernaryOp>(node)->getCondition();

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
	uint32_t MpOperationTernary<T>::optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const
	{
		std::shared_ptr<AstTernaryOp> ternaryNode = std::static_pointer_cast<AstTernaryOp>(node);
		std::shared_ptr<AstNode> branchCond = ternaryNode->getCondition();
		if (branchCond->isImm())
		{
			// optimize an immediate condition
			bool conditionIsTrue = std::static_pointer_cast<AstImm>(branchCond)->getValue<cplx_t>() != cplx_t({ 0, 0 });

			std::shared_ptr<AstNode> nodeOptimized;

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

			nodeOptimized->_parent.reset();
			ternaryNode->getParent()->replaceNode(ternaryNode, nodeOptimized);
		}

		return ErrorCode::kErrorOk;

	}

	// Variable declaration
	template<typename T>
	class MpOperationVarDeclaration : public MpOperation
	{
	public:
		MpOperationVarDeclaration() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};

	template<>
	MpOperationVarDeclaration<double>::MpOperationVarDeclaration() noexcept
		: MpOperation(Signature(1, Signature::type::real), MpOperation::RighttoLeft | MpOperation::IsAssignment, 15)
	{
	}
	template<>
	MpOperationVarDeclaration<cplx_t>::MpOperationVarDeclaration() noexcept
		: MpOperation(Signature(1, Signature::type::complex), MpOperation::RighttoLeft | MpOperation::IsAssignment, 15)
	{
	}

	template<typename T>
	JitVar MpOperationVarDeclaration<T>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		JitVar result;
		std::shared_ptr<AstVarDecl> varDecl = std::static_pointer_cast<AstVarDecl>(node);

		if (varDecl->hasChild())
			result = jc->onNode(varDecl->getChild());

		std::shared_ptr<AstSymbol> sym = varDecl->getSymbol();
		uint32_t slotId = sym->getVarSlotId();

		result.setRO();
		jc->varSlots[slotId] = result;

		return result;
	}

	template<typename T>
	uint32_t MpOperationVarDeclaration<T>::optimize(AstOptimizer * opt, std::shared_ptr<AstNode> node) const
	{
		std::shared_ptr<AstVarDecl> varDecl;
		if (node->getNodeType() == AstNodeType::kAstNodeVarDecl)
		{
			varDecl = std::static_pointer_cast<AstVarDecl>(node);
		}
		else
		{
			return ErrorCode::kErrorInvalidState;
		}

		std::shared_ptr<AstSymbol> sym = varDecl->getSymbol();


		if (varDecl->hasChild())
		{
			std::shared_ptr<AstNode> child = varDecl->getChild();

			if (child->returnsComplex())
			{
				sym->setSymbolFlag(AstSymbolFlags::kAstSymbolIsComplex);
			}
			if (child->isImm())
			{
				sym->setValue(std::static_pointer_cast<AstImm>(child)->getValue<T>());
				sym->setAssigned();
			}
		}

		return ErrorCode::kErrorOk;
	}

	// Assignment
	template<typename T>
	class MpOperationAssignment : public MpOperation
	{
	public:
		MpOperationAssignment() noexcept;
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;
	};

	template<>
	MpOperationAssignment<double>::MpOperationAssignment() noexcept
		: MpOperation(Signature(2, Signature::type::real), MpOperation::RighttoLeft | MpOperation::IsAssignment, 15)
	{
	}
	template<>
	MpOperationAssignment<cplx_t>::MpOperationAssignment()  noexcept
		: MpOperation(Signature(2, Signature::type::complex), MpOperation::RighttoLeft | MpOperation::IsAssignment, 15)
	{
	}

	template<typename T>
	JitVar MpOperationAssignment<T>::compile(JitCompiler * jc, std::shared_ptr<AstNode> node) const
	{
		std::shared_ptr<AstVar> varNode = std::static_pointer_cast<AstVar>(node->getAt(0));
		MATHPRESSO_ASSERT(varNode->getNodeType() == AstNodeType::kAstNodeVar);

		std::shared_ptr<AstSymbol> sym = varNode->getSymbol();
		uint32_t slotId = sym->getVarSlotId();

		JitVar result = jc->onNode(node->getAt(1));
		result.setRO();

		sym->setAltered();
		jc->varSlots[slotId] = result;

		return result;
	}

	template<typename T>
	uint32_t MpOperationAssignment<T>::optimize(AstOptimizer * opt, std::shared_ptr<AstNode> node) const
	{
		std::shared_ptr<AstSymbol> sym = std::static_pointer_cast<AstVar>(node->getAt(0))->getSymbol();
		if (node->getAt(1)->isImm())
		{
			if (sym->isAssigned())
			{
				sym->setValue(std::static_pointer_cast<AstImm>(node->getAt(1))->getValue<T>());
			}
		}
		sym->setAltered();
		return ErrorCode::kErrorOk;
	}

	uint32_t addBuiltinMpObjects(Context * ctx)
	{
		ctx->addObject("+", std::make_shared<MpOperationAdd<double>>());
		ctx->addObject("+", std::make_shared<MpOperationAdd<cplx_t>>());
		ctx->addObject("-", std::make_shared<MpOperationSub<double>>());
		ctx->addObject("-", std::make_shared<MpOperationSub<cplx_t>>());
		ctx->addObject("*", std::make_shared<MpOperationMul<double>>());
		ctx->addObject("*", std::make_shared<MpOperationMul<cplx_t>>());
		ctx->addObject("/", std::make_shared<MpOperationDiv<double>>());
		ctx->addObject("/", std::make_shared<MpOperationDiv<cplx_t>>());
		ctx->addObject("==", std::make_shared<MpOperationEq<double>>());
		ctx->addObject("==", std::make_shared<MpOperationEq<cplx_t>>());
		ctx->addObject("!=", std::make_shared<MpOperationNe<double>>());
		ctx->addObject("!=", std::make_shared<MpOperationNe<cplx_t>>());
		ctx->addObject(">=", std::make_shared<MpOperationGe>());
		ctx->addObject(">", std::make_shared<MpOperationGt>());
		ctx->addObject("<=", std::make_shared<MpOperationLe>());
		ctx->addObject("<", std::make_shared<MpOperationLt>());
		ctx->addObject("min", std::make_shared<MpOperationMin>());
		ctx->addObject("max", std::make_shared<MpOperationMax>());
		ctx->addObject("%", std::make_shared<MpOperationModulo>());
		ctx->addObject("_ternary_", std::make_shared<MpOperationTernary<double>>());
		ctx->addObject("_ternary_", std::make_shared<MpOperationTernary<cplx_t>>());
		ctx->addObject("=", std::make_shared<MpOperationVarDeclaration<double>>());
		ctx->addObject("=", std::make_shared<MpOperationVarDeclaration<cplx_t>>());
		ctx->addObject("=", std::make_shared<MpOperationAssignment<double>>());
		ctx->addObject("=", std::make_shared<MpOperationAssignment<cplx_t>>());
		ctx->addObject("isfinite", std::make_shared<MpOperationIsFinite<double>>());
		ctx->addObject("isfinite", std::make_shared<MpOperationIsFinite<cplx_t>>());
		ctx->addObject("isinf", std::make_shared<MpOperationIsInfinite<double>>());
		ctx->addObject("isinf", std::make_shared<MpOperationIsInfinite<cplx_t>>());
		ctx->addObject("isnan", std::make_shared<MpOperationIsNan<double>>());
		ctx->addObject("isnan", std::make_shared<MpOperationIsNan<cplx_t>>());

		ctx->addObject("-", std::make_shared<MpOperationNeg<double>>());
		ctx->addObject("-", std::make_shared<MpOperationNeg<cplx_t>>());
		ctx->addObject("abs", std::make_shared<MpOperationAbs>());
		_ADD_FUNRC("abs", std::abs);

		ctx->addObject("avg", std::make_shared<MpOperationAvg<double>>());
		ctx->addObject("avg", std::make_shared<MpOperationAvg<cplx_t>>());
		ctx->addObject("recip", std::make_shared<MpOperationRecip<double>>());
		ctx->addObject("recip", std::make_shared<MpOperationRecip<cplx_t>>());
		ctx->addObject("!", std::make_shared<MpOperationNot<double>>());
		ctx->addObject("!", std::make_shared<MpOperationNot<cplx_t>>());

		ctx->addObject("real", std::make_shared<MpOperationGetReal>());
		ctx->addObject("imag", std::make_shared<MpOperationGetImag>());
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

		_ADD_FUNC1("sin", std::sin);
		_ADD_FUNC1("cos", std::cos);
		_ADD_FUNC1("tan", std::tan);
		_ADD_FUNC1("sinh", std::sinh);
		_ADD_FUNC1("cosh", std::cosh);
		_ADD_FUNC1("tanh", std::tanh);
		_ADD_FUNC1("acos", std::acos);
		_ADD_FUNC1("asin", std::asin);
		_ADD_FUNC1("atan", std::atan);
		_ADD_FUNC1("log", std::log);
		_ADD_FUNC1("log10", std::log10);
		_ADD_FUNC1("exp", std::exp);
		_ADD_FUNC1("log2", std::log2);
		_ADD_FUNC2("pow", std::pow);
		ctx->addObject("atan2", _OBJ(static_cast<double(*)(double, double)>(std::atan2)));
		ctx->addObject("hypot", _OBJ(static_cast<double(*)(double, double)>(std::hypot)));
		ctx->addObject("sqrtc", _OBJ(&sqrtc));
		ctx->addObject("sqrtc", _OBJ(static_cast<cplx_t(*)(const cplx_t&)>(std::sqrt)));

		ctx->addObject("_none_", std::make_shared<MpOperationFunc<double, double>>(nullptr, 0));
		ctx->addObject("_none_", std::make_shared<MpOperationFunc<cplx_t, cplx_t>>(nullptr, 0));

		ctx->addObject("?", std::make_shared<MpOperationFunc<double, double>>(nullptr, 2, MpOperation::RighttoLeft, 15));
		ctx->addObject("?", std::make_shared<MpOperationFunc<cplx_t, cplx_t>>(nullptr, 2, MpOperation::RighttoLeft, 15));
		ctx->addObject(":", std::make_shared<MpOperationFunc<double, double>>(nullptr, 2, MpOperation::RighttoLeft, 15));
		ctx->addObject(":", std::make_shared<MpOperationFunc<cplx_t, cplx_t>>(nullptr, 2, MpOperation::RighttoLeft, 15));

		ctx->addConstant("NaN", mpGetNan());
		ctx->addConstant("INF", mpGetInf());
		ctx->addConstant("PI", 3.14159265358979323846);
		ctx->addConstant("E", 2.7182818284590452354);
		ctx->addConstant("i", { 0, 1 });

		return 0;
	}

} // end namespace mathpresso
