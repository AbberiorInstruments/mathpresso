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

	// MpOperationFunc
	JitVar MpOperationFunc::compile(JitCompiler* jc, AstNode * node) 
	{
		asmjit::X86Xmm result = node->returnsComplex() ? jc->cc->newXmmPd() : jc->cc->newXmmSd();
		asmjit::X86Xmm args[8];

		if (!node->takesComplex())
		{
			bool returnsComplex = (flags_ & OpFlagDReturnsC) != 0;
			if (node->returnsComplex() != returnsComplex  || !fnD_)
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
			bool returnsComplex = (flags_ & OpFlagCReturnsD) == 0;
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
				std::vector<std::complex<double>> args;
				for (size_t i = 0; i < nargs_; i++)
				{
					args.push_back((static_cast<AstImm*>(node->getAt(i)))->getValueCplx());
				}

				if (b_returns_complex)
				{
					ret->setValue(((mpFuncpCtoC)fnC_)(args.data()));
				}
				else
				{
					ret->setValue(((mpFuncpCtoD)fnC_)(args.data()));
				}
			}
			else
			{
				std::vector<double> args;
				for (size_t i = 0; i < nargs_; i++)
				{
					args.push_back((static_cast<AstImm*>(node->getAt(i)))->getValue());
				}

				if (b_returns_complex)
				{
					ret->setValue(((mpFuncpDtoC)fnD_)(args.data()));
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

		return ErrorCode::kErrorOk;
	}

	double MpOperationFunc::fnDouble(AstOptimizer * opt, AstNode * node) {
		return 0.0;
	}

	std::complex<double> MpOperationFunc::fnComplex(AstOptimizer * opt, AstNode * node) {
		return std::complex<double>();
	}

	void MpOperationFunc::setFn(void* fn, bool isComplex) {
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
	
	// MpOperationFuncAsm
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

				if (asmD_)
				{
					std::vector<JitVar> args;

					for (size_t i = 0; i < node->getLength(); i++)
					{
						args.push_back(jc->registerVar(jc->onNode(node->getAt(i))));
					}

					return asmD_(jc, args.data());
				}
				else
				{// if by any chance no assembler-generator is available:
					asmjit::X86Xmm result = node->returnsComplex() ? jc->cc->newXmmPd() : jc->cc->newXmmSd();
					asmjit::X86Xmm args[8];

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

					return JitVar(result, JitVar::FLAG_NONE);
				}
				
			}
			else
			{
				bool returnsComplex = (flags_ & OpFlagCReturnsD) == 0;
				if (node->returnsComplex() != returnsComplex || !fnC_)
				{
					// Should never happen, as the optimizer should have taken care of that. Remove later
					throw std::runtime_error("Implementation error!");
				}

				if (asmC_)
				{
					std::vector<JitVar> args;
					for (size_t i = 0; i < node->getLength(); i++)
					{
						args.push_back(jc->registerVarComplex(jc->onNode(node->getAt(i))));
					}

					return asmC_(jc, args.data());
				}
				else
				{ // if by any chance no assembler-generator is available:
					asmjit::X86Xmm result = node->returnsComplex() ? jc->cc->newXmmPd() : jc->cc->newXmmSd();
					asmjit::X86Xmm args[8];
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
					return JitVar(result, JitVar::FLAG_NONE);
				}
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

	// mpOperationBinary
	JitVar MpOperationBinary::compile(JitCompiler* jc, AstNode * node) 
	{
		JitVar vl, vr;
		AstNode* left = static_cast<AstBinaryOp*>(node)->getLeft();
		AstNode* right = static_cast<AstBinaryOp*>(node)->getRight();

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
		if ((needs_complex || hasFlag(OpHasNoReal)) && !hasFlag(OpHasNoComplex))
		{
			node->addNodeFlags(AstNodeFlags::kAstReturnsComplex | AstNodeFlags::kAstTakesComplex);
		}

		if (lIsImm && rIsImm)
		{
			// optimize a calculation with two immediates.
			AstImm* lNode = static_cast<AstImm*>(left);
			AstImm* rNode = static_cast<AstImm*>(right);
			if (needs_complex && !hasFlag(OpHasNoComplex))
			{
				lNode->setValue(calculateComplex(lNode->getValueCplx(), rNode->getValueCplx()));
			}
			else if (!needs_complex && !hasFlag(OpHasNoReal))
			{
				lNode->setValue(calculateReal(lNode->getValue(), rNode->getValue()));
			}
			else if (!needs_complex && !hasFlag(OpHasNoComplex))
			{
				lNode->setValue(calculateComplex(lNode->getValueCplx(), rNode->getValueCplx()));
			}
			else
			{
				throw std::runtime_error("Wrong implementation of MpOperationBinary!");
			}
			
			// setValue sets the correct flags automatically.
			static_cast<AstBinaryOp*>(node)->unlinkLeft();
			node->getParent()->replaceNode(node, lNode);
			opt->getAst()->deleteNode(node);
		}
		else if (lIsImm)
		{
			AstImm* lNode = static_cast<AstImm*>(left);
			// if the node is real, the imaginary part is set to zero by default.
			if ((hasFlag(OpFlagNopIfLZero) && lNode->getValueCplx() == std::complex<double>(0.0, 0.0)) ||
				(hasFlag(OpFlagNopIfLOne) && lNode->getValueCplx() == std::complex<double>(1.0, 0.0)))
			{
				static_cast<AstBinaryOp*>(node)->unlinkRight();
				node->getParent()->replaceNode(node, right);
				opt->getAst()->deleteNode(node);
			}
		}
		else if (rIsImm)
		{
			AstImm* rNode = static_cast<AstImm*>(right);
			
			if ((hasFlag(OpFlagNopIfRZero) && rNode->getValueCplx() == std::complex<double>(0.0, 0.0)) ||
				(hasFlag(OpFlagNopIfROne) &&  rNode->getValueCplx() == std::complex<double>(1.0, 0.0)))
			{
				static_cast<AstBinaryOp*>(node)->unlinkLeft();
				node->getParent()->replaceNode(node, left);
				opt->getAst()->deleteNode(node);
			}
			
		}
		return ErrorCode::kErrorOk;
	}


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

	double MpOperationAdd::calculateReal(double vl, double vr) {	return vl + vr;	}
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
		jc->cc->shufpd(vr.getXmm(), vr.getXmm(), 1);
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
		jc->cc->shufpd(vr.getXmm(), vr.getXmm(), 1);
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

	double MpOperationNe::calculateReal(double vl, double vr) { return vl != vr ? 1.0: 0.0; }
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

	JitVar MpOperationLt::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		throw std::runtime_error("Wrong implementation of MpOperationbinary!");
	}

	double MpOperationLt::calculateReal(double vl, double vr) { return vl < vr ? 1.0 : 0.0; }
	std::complex<double> MpOperationLt::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return std::complex<double>(NAN, NAN); }

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

	JitVar MpOperationLe::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		throw std::runtime_error("Wrong implementation of MpOperationbinary!");
	}

	double MpOperationLe::calculateReal(double vl, double vr) { return vl <= vr ? 1.0 : 0.0; }
	std::complex<double> MpOperationLe::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return std::complex<double>(NAN, NAN); }

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

	JitVar MpOperationGt::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		throw std::runtime_error("Wrong implementation of MpOperationbinary!");
	}

	double MpOperationGt::calculateReal(double vl, double vr) { return vl > vr ? 1.0 : 0.0; }
	std::complex<double> MpOperationGt::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return std::complex<double>(NAN, NAN); }

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

	JitVar MpOperationGe::generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr)
	{
		throw std::runtime_error("Wrong implementation of MpOperationbinary!");
	}

	double MpOperationGe::calculateReal(double vl, double vr) { return vl >= vr ? 1.0 : 0.0; }
	std::complex<double> MpOperationGe::calculateComplex(std::complex<double> vl, std::complex<double> vr) { return std::complex<double>(NAN, NAN); }

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
		if (static_cast<AstBinaryOp*>(node)->getRight() != lastColon) {
			// remove left branch from the AST.
			branchLeft = static_cast<AstBinaryOp*>(node)->getRight();
			static_cast<AstBinaryOp*>(node)->setRight(nullptr);
			branchLeft->_parent = nullptr;

			// correct the right path.
			AstBinaryOp* preLastColon = static_cast<AstBinaryOp*>(lastColon->getParent());
			preLastColon->replaceAt(1, lastColon->getLeft());

		}
		// i.e.: cond1 ? a : b
		else {
			// remove left branch from the AST.
			lastColon->setLeft(nullptr);
			branchLeft->_parent = nullptr;
		}

		// create the new Ternary Node.
		AstTernaryOp* ternaryNode = node->getAst()->newNode<AstTernaryOp>(kOpQMark);
		ternaryNode->setCondition(branchCondition);
		ternaryNode->setLeft(branchLeft);
		ternaryNode->setRight(branchRight);

		AstBinaryOp* oldNode = static_cast<AstBinaryOp*>(node);

		// add the new node to the AST.
		node->getParent()->replaceNode(node, ternaryNode);

		// clean up:
		lastColon->setLeft(nullptr);
		opt->_ast->deleteNode(lastColon);
		opt->_ast->deleteNode(node);

		MATHPRESSO_PROPAGATE(opt->onNode(ternaryNode->getCondition()));
		AstNode* branchCond = ternaryNode->getCondition();
		if (branchCond->isImm()) {
			// optimize an immediate condition
			bool conditionIsTrue = static_cast<AstImm*>(branchCond)->getValueCplx() != std::complex<double>({ 0, 0 });

			AstNode* nodeOptimized;

			if (conditionIsTrue) {
				nodeOptimized = ternaryNode->getLeft();
				ternaryNode->setLeft(nullptr);
			}
			else {
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
		else {
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
	JitVar mathpresso::MpOperationAssignment::compile(JitCompiler * jc, AstNode * node) {
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

		if (varDecl->hasChild()) {
			MATHPRESSO_PROPAGATE(opt->onNode(varDecl->getChild()));
			AstNode* child = varDecl->getChild();

			if (child->isImm())
			{
				if (child->returnsComplex())
				{
					sym->setValue(static_cast<AstImm*>(child)->getValueCplx());
					varDecl->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);
					sym->setSymbolFlag(kAstSymbolIsComplex);
				}
				else {
					sym->setValue(static_cast<AstImm*>(child)->getValue());
				}
				sym->setAssigned();
			}
		}

		return kErrorOk;
	}


} // end namespace mathpresso
