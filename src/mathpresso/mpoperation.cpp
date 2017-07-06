// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

#include "mpoperation_p.h"
#include "mpast_p.h"
#include "mpcompiler_p.h"
#include "mpoptimizer_p.h"
#include "asmjit\x86\x86operand.h"
#include "asmjit\x86\x86inst.h"

#include <complex>

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
				return ErrorCode::kErrorInvalidArgument;

			node->addNodeFlags(AstNodeFlags::kAstTakesComplex);
			b_returns_complex = (flags_ & OpFlagCReturnsD) == 0;
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

		// optimize all-immediate calls:
		if (b_all_imm) 
		{
			AstImm* ret = opt->getAst()->newNode<AstImm>(0);
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
			node->getParent()->replaceNode(node, ret);
			opt->onNode(ret);

			opt->getAst()->deleteNode(node);
		}

		return ErrorCode::kErrorOk;
	}

	void MpOperationFunc::setFn(void* fn, bool isComplex) {
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

	// -- MpOperationFuncAsm

	JitVar mathpresso::MpOperationFuncAsm::compile(JitCompiler * jc, AstNode * node) {
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
				bool returnsComplex = (flags_ & OpFlagCReturnsD) == 0;
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

	//-- MpOperationAdd

	JitVar MpOperationBinary::compile(JitCompiler* jc, AstNode * node) 
	{
		JitVar vl, vr;
		AstNode* left = static_cast<AstBinaryOp*>(node)->getLeft();
		AstNode* right = static_cast<AstBinaryOp*>(node)->getRight();

		if (node->takesComplex())
		{
			// check whether the vars are the same, to reduce memory-operations
			if (left->isVar() && right->isVar() && hasFlag(OpIsCommutativ) &&
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
			return compComplex(jc, vl, vr);
		}
		else
		{
			// check whether the vars are the same, to reduce memory-operations
			if (left->isVar() && right->isVar() && hasFlag(OpIsCommutativ) &&
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
			return compReal(jc, vl, vr);
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
		if (needs_complex)
		{
			node->addNodeFlags(AstNodeFlags::kAstReturnsComplex | AstNodeFlags::kAstTakesComplex);
		}

		if (lIsImm && rIsImm)
		{
			// clear an immediate addition.
			AstImm* lNode = static_cast<AstImm*>(left);
			AstImm* rNode = static_cast<AstImm*>(right);
			if (needs_complex)
			{
				lNode->setValue(optComplex(lNode->getValueCplx(), rNode->getValueCplx()));
			}
			else
			{
				lNode->setValue(optReal(lNode->getValue(), rNode->getValue()));
			}
			// setValue sets the correct flags automaticaly.

			static_cast<AstBinaryOp*>(node)->unlinkLeft();
			node->getParent()->replaceNode(node, lNode);
			opt->getAst()->deleteNode(node);
		}
		else if (lIsImm)
		{
			AstImm* lNode = static_cast<AstImm*>(left);
			// if the node is real, the imaginary part is set to zero by default.
			if ((lNode->getValueCplx() == std::complex<double>(0.0, 0.0) && hasFlag(OpFlagNopIfLZero)) ||
				(lNode->getValueCplx() == std::complex<double>(1.0, 0.0) && hasFlag(OpFlagNopIfLOne)))
			{
				static_cast<AstBinaryOp*>(node)->unlinkRight();
				node->getParent()->replaceNode(node, right);
				opt->getAst()->deleteNode(node);
			}
		}
		else if (rIsImm) {
			AstImm* rNode = static_cast<AstImm*>(right);
			
			if ((rNode->getValueCplx() == std::complex<double>(0.0, 0.0) && hasFlag(OpFlagNopIfRZero)) ||
				(rNode->getValueCplx() == std::complex<double>(1.0, 0.0) && hasFlag(OpFlagNopIfROne)))
			{
				static_cast<AstBinaryOp*>(node)->unlinkLeft();
				node->getParent()->replaceNode(node, left);
				opt->getAst()->deleteNode(node);
			}
			
		}
		return ErrorCode::kErrorOk;
	}


	// Addition
	JitVar MpOperationAdd::compReal(JitCompiler * jc, JitVar vl, JitVar vr) 
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

	JitVar MpOperationAdd::compComplex(JitCompiler * jc, JitVar vl, JitVar vr) 
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

	double MpOperationAdd::optReal(double vl, double vr) {	return vl + vr;	}

	std::complex<double> MpOperationAdd::optComplex(std::complex<double> vl, std::complex<double> vr) { return vl + vr; }


	//Substraction
	JitVar MpOperationSub::compReal(JitCompiler * jc, JitVar vl, JitVar vr) {
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

	JitVar MpOperationSub::compComplex(JitCompiler * jc, JitVar vl, JitVar vr) {
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

	double MpOperationSub::optReal(double vl, double vr) { return vl - vr; }

	std::complex<double> MpOperationSub::optComplex(std::complex<double> vl, std::complex<double> vr) { return vl - vr; }

	// Multiplication
	JitVar MpOperationMul::compReal(JitCompiler * jc, JitVar vl, JitVar vr) {
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

	JitVar MpOperationMul::compComplex(JitCompiler * jc, JitVar vl, JitVar vr) {
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
		jc->cc->shufpd(vr.getXmm(), vr.getXmm(), 1);
		jc->cc->pxor(vr.getXmm(), negateImag.getMem());
		jc->cc->mulpd(vl.getXmm(), vr.getXmm());
		jc->cc->hsubpd(ret.getXmm(), vl.getXmm());
		return ret;
	}

	double MpOperationMul::optReal(double vl, double vr) { return vl * vr; }

	std::complex<double> MpOperationMul::optComplex(std::complex<double> vl, std::complex<double> vr) { return vl * vr; }

	// Division
	JitVar MpOperationDiv::compReal(JitCompiler * jc, JitVar vl, JitVar vr) {
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

	JitVar MpOperationDiv::compComplex(JitCompiler * jc, JitVar vl, JitVar vr) {
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
		jc->cc->shufpd(vr.getXmm(), vr.getXmm(), 1);
		jc->cc->pxor(vr.getXmm(), negateImag.getMem());
		jc->cc->mulpd(vl.getXmm(), vr.getXmm());
		jc->cc->hsubpd(ret.getXmm(), vl.getXmm());

		jc->cc->mulpd(vr.getXmm(), vr.getXmm());
		jc->cc->haddpd(vr.getXmm(), vr.getXmm());
		jc->cc->divpd(ret.getXmm(), vr.getXmm());
		return ret;
	}

	double MpOperationDiv::optReal(double vl, double vr) { return vl / vr; }

	std::complex<double> MpOperationDiv::optComplex(std::complex<double> vl, std::complex<double> vr) { return vl / vr; }

} // end namespace mathpresso
