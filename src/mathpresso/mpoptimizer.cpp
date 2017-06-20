// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include "./mpast_p.h"
#include "./mpeval_p.h"
#include "./mpoptimizer_p.h"

namespace mathpresso {

	// ============================================================================
	// [mpsl::AstOptimizer - Construction / Destruction]
	// ============================================================================

	AstOptimizer::AstOptimizer(AstBuilder* ast, ErrorReporter* errorReporter)
		: AstVisitor(ast),
		_errorReporter(errorReporter) {}
	AstOptimizer::~AstOptimizer() {}

	// ============================================================================
	// [mpsl::AstOptimizer - OnNode]
	// ============================================================================

	Error AstOptimizer::onBlock(AstBlock* node) {
		// Prevent removing nodes that are not stored in pure `AstBlock`. For example
		// function call inherits from `AstBlock`, but it needs each expression passed.
		bool alterable = node->getNodeType() == kAstNodeBlock;

		uint32_t i = 0;
		uint32_t curCount = node->getLength();
		uint32_t oldCount;
		bool isComplex = false;

		while (i < curCount) 
		{
			MATHPRESSO_PROPAGATE(onNode(node->getAt(i)));

			oldCount = curCount;
			curCount = node->getLength();

			if (curCount < oldCount) {
				if (!alterable)
					return MATHPRESSO_TRACE_ERROR(kErrorInvalidState);
				continue;
			}

			if (alterable && (node->getAt(i)->isImm())) 
			{
				_ast->deleteNode(node->removeAt(i));
				curCount--;
				continue;
			}
			// TODO: check semantics here
			isComplex |= node->getAt(i)->returnsComplex();
			i++;
		}

		if (isComplex)
		{
			node->addNodeFlags(kAstReturnsComplex | kAstTakesComplex);
		}

		return kErrorOk;
	}

	Error AstOptimizer::onVarDecl(AstVarDecl* node) {
		AstSymbol* sym = node->getSymbol();

		if (node->hasChild()) {
			MATHPRESSO_PROPAGATE(onNode(node->getChild()));
			AstNode* child = node->getChild();

			if (child->isImm())
			{
				if (child ->returnsComplex())
					sym->setValue(static_cast<AstImm*>(child)->getValue());
				else
					sym->setValue(static_cast<AstImm*>(child)->getValueComp());
			}
		}

		return kErrorOk;
	}

	Error AstOptimizer::onVar(AstVar* node) {
		AstSymbol* sym = node->getSymbol();
		bool b_complex = node ->returnsComplex();

		if (sym->isAssigned() && !node->hasNodeFlag(kAstNodeHasSideEffect)) 
		{
			AstImm* imm;

			if (!b_complex)
			{
				imm = _ast->newNode<AstImm>(sym->getValue());
			}
			else
			{
				imm = _ast->newNode<AstImm>(sym->getValueComp());
			}
			_ast->deleteNode(node->getParent()->replaceNode(node, imm));
		}
		return kErrorOk;
	}

	Error AstOptimizer::onImm(AstImm* node) {
		return kErrorOk;
	}

	Error AstOptimizer::onUnaryOp(AstUnaryOp* node)
	{
		const OpInfo& op = OpInfo::get(node->getOp());

		if (op.isComplex())
		{
			node->addNodeFlags(kAstTakesComplex);
		}

		if (op.returnsComplex())
		{	
			node->addNodeFlags(kAstReturnsComplex);
		}

		MATHPRESSO_PROPAGATE(onNode(node->getChild()));
		AstNode* child = node->getChild();
		//if (child->hasNodeFlag(kAstComplex) && !node->hasNodeFlag(kAstComplex)) {
	  //	  return _errorReporter->onError(kErrorInvalidArgument, node->getPosition(),
	  //		  "Operator '%s' wants a non-complex Parameter, but gets a complex one.", op.name);
		//}

		if (child->isImm())
		{
			AstImm* p_imm = static_cast<AstImm*>(node->getChild());

			if (!p_imm ->returnsComplex())
			{
				if (!op.returnsComplex())
				{
					double value = p_imm->getValue();

					if (!op.isIntrinsic()) {
						switch (node->getOp()) {
						case kOpNeg: value = -value; break;
						case kOpNot: value = (value == 0); break;

						default:
							return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
								"Invalid unary operation '%s'.", op.name);
						}
					}
					else {
						if (op.funcD)
							value = ((Arg1Func)op.funcD)(value);
						else {
							return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
								"Invalid unary operation '%s'.", op.name);
						}
						/*
						case kOpIsNan: value = mpIsNan(value); break;
						case kOpIsInf: value = mpIsInf(value); break;
						case kOpIsFinite: value = mpIsFinite(value); break;
						case kOpSignBit: value = mpSignBit(value); break;

						case kOpRound: value = mpRound(value); break;
						case kOpRoundEven: value = mpRoundEven(value); break;
						case kOpTrunc: value = mpTrunc(value); break;
						case kOpFloor: value = mpFloor(value); break;
						case kOpCeil: value = mpCeil(value); break;

						case kOpAbs: value = mpAbs(value); break;
						case kOpExp: value = mpExp(value); break;

						case kOpLog: value = mpLog(value); break;
						case kOpLog2: value = mpLog2(value); break;
						case kOpLog10: value = mpLog10(value); break;

						case kOpSqrt: value = mpSqrt(value); break;
						case kOpFrac: value = mpFrac(value); break;
						case kOpRecip: value = mpRecip(value); break;

						case kOpSin: value = mpSin(value); break;
						case kOpCos: value = mpCos(value); break;
						case kOpTan: value = mpTan(value); break;

						case kOpSinh: value = mpSinh(value); break;
						case kOpCosh: value = mpCosh(value); break;
						case kOpTanh: value = mpTanh(value); break;

						case kOpAsin: value = mpAsin(value); break;
						case kOpAcos: value = mpAcos(value); break;
						case kOpAtan: value = mpAtan(value); break;
						*/
					}

					p_imm->setValue(value);

					node->unlinkChild();
					node->getParent()->replaceNode(node, p_imm);
				}
				else
				{
					double operand = static_cast<AstImm*>(node->getChild())->getValue();
					std::complex<double> value(0, 0);
					if (op.funcC) {
						value = ((mpFuncpCtoC)op.funcC)(&std::complex<double>(operand, 0));
					} else {
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"invalid unary operation '%s'.", op.name);
					}

					AstImm* replacement = _ast->newNode<AstImm>(value);
					node->getParent()->replaceNode(node, replacement);
				}
				_ast->deleteNode(node);
			}
			else
			{
				std::complex<double> value = p_imm->getValueComp();

				if (op.funcC) 
				{
					if ((op.flags & kOpFlagComplex) != 0) 
					{
						p_imm->setValue(((mpFuncpCtoD)op.funcC)(&value));
						p_imm->removeNodeFlags(kAstReturnsComplex);
						//case kOpReal: value = mpGetReal(&value); p_imm->removeNodeFlags(kAstReturnsComplex);
						//case kOpImag: value = mpGetImag(&value); p_imm->removeNodeFlags(kAstReturnsComplex);
					}
					else 
					{
						p_imm->setValue(((mpFuncpCtoC)op.funcC)(&value));
						/*
						case kOpExp: value = mpExpC(&value); break;

						case kOpLog: value = mpLogC(&value); break;
						case kOpLog2: value = mpLog2C(&value); break;
						case kOpLog10: value = mpLog10C(&value); break;

						case kOpSin: value = mpSinC(&value); break;
						case kOpCos: value = mpCosC(&value); break;
						case kOpTan: value = mpTanC(&value); break;

						case kOpSinh: value = mpSinhC(&value); break;
						case kOpCosh: value = mpCoshC(&value); break;
						case kOpTanh: value = mpTanhC(&value); break;

						case kOpAsin: value = mpAsinC(&value); break;
						case kOpAcos: value = mpAcosC(&value); break;
						case kOpAtan: value = mpAtanC(&value); break;
						*/
					}
				}
				else 
				{
					return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
						"Invalid unary operation with complex parameters '%s'.", op.name);
				}
				
				if (p_imm->returnsComplex()) {
					node->unlinkChild();
					node->getParent()->replaceNode(node, p_imm);
					_ast->deleteNode(node);
					onNode(p_imm);
				}
				else {
					AstImm* newChild = _ast->newNode<AstImm>(value.real());
					node->getParent()->replaceNode(node, newChild);
					_ast->deleteNode(node);
					onNode(newChild);
				}
			}
		}
		else if (child->getNodeType() == kAstNodeUnaryOp && node->getOp() == child->getOp())
		{
			// Simplify `-(-(x))` -> `x`.
			if (node->getOp() == kOpNeg || node->getOp() == kOpConjug)
			{
				AstNode* childOfChild = static_cast<AstUnaryOp*>(child)->unlinkChild();
				node->getParent()->replaceNode(node, childOfChild);
				_ast->deleteNode(node);
			}
		}
		else if (node->getChild()->returnsComplex()) 
		{
			// If it is not found, also check for a user function with the same name and 
			// arguments
			switch (node->getOp()) {
			case kOpLog: node->setOp(kOpLogC); break;
			case kOpLog2: node->setOp(kOpLog2C); break;
			case kOpLog10: node->setOp(kOpLog10C); break;
			case kOpSin: node->setOp(kOpSinC); break;
			case kOpCos: node->setOp(kOpCosC); break;
			case kOpTan: node->setOp(kOpTanC); break;
			case kOpSinh: node->setOp(kOpSinhC); break;
			case kOpCosh: node->setOp(kOpCoshC); break;
			case kOpTanh: node->setOp(kOpTanhC); break;
			case kOpAsin: node->setOp(kOpAsinC); break;
			case kOpAcos: node->setOp(kOpAcosC); break;
			case kOpAtan: node->setOp(kOpAtanC); break;
			case kOpExp: node->setOp(kOpExpC); break;
			default: return kErrorOk;
			}
			
			node->addNodeFlags(kAstReturnsComplex | kAstTakesComplex);
		}

		return kErrorOk;
	}

	Error AstOptimizer::onBinaryOp(AstBinaryOp* node)
	{
		const OpInfo& op = OpInfo::get(node->getOp());

		AstNode* left = node->getLeft();
		AstNode* right = node->getRight();

		if (op.isAssignment())
			left->addNodeFlags(kAstNodeHasSideEffect);

		if (node->getOp() == kOpColon) {
			return _errorReporter->onError(kErrorInvalidSyntax, node->getPosition(),
				"Invalid Operator ':' found. No '?' to complete the ternary expression.");
		}

		if (node->getOp() == kOpQMark) {

			AstBinaryOp* lastColon = node;
			// go to the last Colon after question-marks.
			while (lastColon->getOp() == kOpQMark) {
				lastColon = static_cast<AstBinaryOp*>(lastColon->getRight());
			}
			while (lastColon->getRight()->getOp() == kOpColon) {
				lastColon = static_cast<AstBinaryOp*>(lastColon->getRight());
			}

			if (lastColon->getOp() != kOpColon) {
				return _errorReporter->onError(kErrorInvalidSyntax, node->getPosition(),
					"Invalid ternary operation. Expected a ':', found '%s' instead.", OpInfo::get(lastColon->getOp()).name);
			}


			AstNode* branchCondition = node->getLeft();
			AstNode* branchLeft = lastColon->getLeft();
			AstNode* branchRight = lastColon->getRight();

			if (!OpInfo::get(branchCondition->getOp()).isCondition()) {
				return _errorReporter->onError(kErrorInvalidArgument, node->getPosition(),
					"Not a condition.");
			}
			// remove branchCondition from the AST
			node->setLeft(NULL);
			branchCondition->_parent = NULL;

			// remove the right path from the AST.
			lastColon->setRight(NULL);
			branchRight->_parent = NULL;


			// Distinguish between a complex and a non-complex case:
			// i.e.: cond1 ? cond2 ? a : b : c
			if (node->getRight() != lastColon) {
				// remove left branch from the AST.
				branchLeft = node->getRight();
				node->setRight(NULL);
				branchLeft->_parent = NULL;

				// correct the right path.
				AstBinaryOp* preLastColon = static_cast<AstBinaryOp*>(lastColon->getParent());
				preLastColon->replaceAt(1, lastColon->getLeft());

			}
			// i.e.: cond1 ? a : b
			else {
				// remove left branch from the AST.
				lastColon->setLeft(NULL);
				branchLeft->_parent = NULL;
			}

			// create the new Ternary Node.
			AstTernaryOp* newNode = node->getAst()->newNode<AstTernaryOp>(kOpQMark);
			newNode->setCondition(branchCondition);
			newNode->setLeft(branchLeft);
			newNode->setRight(branchRight);

			AstBinaryOp* oldNode = node;

			// add the new node to the AST.
			node->getParent()->replaceNode(node, newNode);

			// clean up:
			lastColon->setLeft(NULL);
			_ast->deleteNode(lastColon);
			_ast->deleteNode(node);


			return onTernaryOp(newNode);
		}

		MATHPRESSO_PROPAGATE(onNode(left));
		left = node->getLeft();

		MATHPRESSO_PROPAGATE(onNode(right));
		right = node->getRight();


		bool lIsImm = left->isImm();
		bool rIsImm = right->isImm();
		bool needs_complex = node->getLeft()->returnsComplex() || node->getRight()->returnsComplex();

		if (!needs_complex)
		{
			if (left->returnsComplex() || right->returnsComplex())
				return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
					"Expect non-complex parameters for operation '%s'.", op.name); // should never be reached

			// If both nodes are values it's easy, just fold them into a single one.
			if (lIsImm && rIsImm)
			{
				AstImm* lNode = static_cast<AstImm*>(left);
				AstImm* rNode = static_cast<AstImm*>(right);

				double lVal = lNode->getValue();
				double rVal = rNode->getValue();
				double result = 0.0;
				if (!op.isIntrinsic())
				{
					switch (node->getOp()) 
					{
					case kOpEq: result = lVal == rVal; break;
					case kOpNe: result = lVal != rVal; break;
					case kOpLt: result = lVal < rVal; break;
					case kOpLe: result = lVal <= rVal; break;
					case kOpGt: result = lVal > rVal; break;
					case kOpGe: result = lVal >= rVal; break;
					case kOpAdd: result = lVal + rVal; break;
					case kOpSub: result = lVal - rVal; break;
					case kOpMul: result = lVal * rVal; break;
					case kOpDiv: result = lVal / rVal; break;
					default:
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid binary operation '%s'.", op.name);
					}
				}
				else
				{
					if (op.funcD) 
					{
						result = ((Arg2Func)op.funcD)(lVal, rVal);
					}
					else 
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid binary operation '%s'.", op.name);
					}

					/*case kOpMod: result = mpMod(lVal, rVal); break;
					case kOpAvg: result = mpAvg(lVal, rVal); break;
					case kOpMin: result = mpMin(lVal, rVal); break;
					case kOpMax: result = mpMax(lVal, rVal); break;
					case kOpPow: result = mpPow(lVal, rVal); break;
					case kOpAtan2: result = mpAtan2(lVal, rVal); break;
					case kOpHypot: result = mpHypot(lVal, rVal); break;
					case kOpCopySign: result = mpCopySign(lVal, rVal); break;*/
				}

				
				

				lNode->setValue(result);
				node->unlinkLeft();
				node->getParent()->replaceNode(node, lNode);

				_ast->deleteNode(node);
			}
			// There is still a little optimization opportunity.
			else if (lIsImm)
			{
				AstImm* lNode = static_cast<AstImm*>(left);
				double val = lNode->getValue();

				if ((val == 0.0 && (op.flags & kOpFlagNopIfLZero)) ||
					(val == 1.0 && (op.flags & kOpFlagNopIfLOne)))
				{
					node->unlinkRight();
					node->getParent()->replaceNode(node, right);

					_ast->deleteNode(node);
				}
			}
			else if (rIsImm) {
				AstImm* rNode = static_cast<AstImm*>(right);
				double val = rNode->getValue();

				// Evaluate an assignment.
				if (op.isAssignment() && left->isVar())
				{
					AstSymbol* sym = static_cast<AstVar*>(left)->getSymbol();
					if (op.type == kOpAssign || sym->isAssigned())
					{
						sym->setValue(val);
						sym->setAssigned();
					}
				}
				else
				{
					if ((val == 0.0 && (op.flags & kOpFlagNopIfRZero)) ||
						(val == 1.0 && (op.flags & kOpFlagNopIfROne)))
					{
						node->unlinkLeft();
						node->getParent()->replaceNode(node, left);

						_ast->deleteNode(node);
					}
				}
			}
		}
		else
		{
			node->addNodeFlags(kAstTakesComplex);
			if (op.returnsComplex() || (!op.returnsComplex() && !op.isComplex()))
				node->addNodeFlags(kAstReturnsComplex);

			// if we have to calculate in complex, and one of the operands is an immediate, it should be converted to complex.
			if (left->isImm()&& !left->takesComplex()) {
				left->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);
			}
			if (right->isImm() && !right->takesComplex()) {
				right->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);
			}

			// complex immediates.
			if (left->isImm() && right->isImm() && left ->takesComplex() && right ->takesComplex()) {
				AstImm* lNode = static_cast<AstImm*>(left);
				AstImm* rNode = static_cast<AstImm*>(right);

				std::complex<double> result;
				if (!op.isIntrinsic()) {
					switch (node->getOp()) {
					case kOpAdd: result = lNode->getValueComp() + rNode->getValueComp(); break;
					case kOpSub: result = lNode->getValueComp() - rNode->getValueComp(); break;
					case kOpMul: result = lNode->getValueComp() * rNode->getValueComp(); break;
					case kOpDiv: result = lNode->getValueComp() / rNode->getValueComp(); break;


					default:
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid complex binary operation '%s'.", op.name);
					}
				}
				else {
					if (op.funcC) {
						result = pow(lNode->getValueComp(), rNode->getValueComp());
					}
					else {
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid complex binary operation '%s'.", op.name);
					}
				}

				rNode->setValue(result);
				node->unlinkRight();
				node->getParent()->replaceNode(node, rNode);

				_ast->deleteNode(node);
				onNode(rNode);
				return kErrorOk;
			}

			if (node->getOp() == kOpPow)
			{
				node->setOp(kOpPowC);
			}

		}

		return kErrorOk;
	}

	Error AstOptimizer::onTernaryOp(AstTernaryOp* node) {
		MATHPRESSO_PROPAGATE(onNode(node->getCondition()));
		AstNode* branchCond = node->getCondition();
		if (branchCond->isImm()) {
			bool conditionIsTrue = static_cast<AstImm*>(branchCond)->getValueComp() != std::complex<double>({ 0, 0 });
			
			AstNode* newNode;

			if (conditionIsTrue) {
				newNode = node->getLeft();
				node->setLeft(nullptr);
			}
			else {
				newNode = node->getRight();
				node->setRight(nullptr);
			}

			newNode->_parent = nullptr;
			node->getParent()->replaceNode(node, newNode);
			node->setCondition(nullptr);
			node->setLeft(nullptr);
			node->setRight(nullptr);
			_ast->deleteNode(node);
			MATHPRESSO_PROPAGATE(onNode(newNode));

		}
		else {
			node->removeNodeFlags(kAstTakesComplex);
			MATHPRESSO_PROPAGATE(onNode(node->getLeft()));
			MATHPRESSO_PROPAGATE(onNode(node->getRight()));
			bool needs_complex = node->getLeft()->returnsComplex() | node->getRight()->returnsComplex();
			if (needs_complex)
			{
				node->addNodeFlags(kAstReturnsComplex | kAstTakesComplex);
			}
		}
		return kErrorOk;

	}

	// TODO: functions that have two versions with different parameter types: 
	Error AstOptimizer::onCall(AstCall* node) 
	{
		AstSymbol* sym = node->getSymbol();
		uint32_t i, count = node->getLength();

		bool b_need_cplx = false;
		for (i = 0; i < count; i++)
		{
			MATHPRESSO_PROPAGATE(onNode(node->getAt(i)));
			b_need_cplx |= node->getAt(i)->takesComplex();
		}
		
		// If we have a complex argument or no real function version, use the function
		// version that takes complex arguments
		if (b_need_cplx || !sym ->getFuncPtr()) 
		{
			node->addNodeFlags(kAstTakesComplex);

			// Need function that takes complex arguments here!
			if (!sym ->getFuncPtr(true))
				return kErrorInvalidArgument;

			for (i = 0; i < count; i++) 
			{
				AstNode *tmp = node->getAt(i);
				tmp->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);
			}

			if (!sym->hasSymbolFlag(kAstSymbolComplexFunctionReturnsReal))
			{
				node->addNodeFlags(kAstReturnsComplex);
			}
		}
		else
		{
			if (sym->hasSymbolFlag(kAstSymbolRealFunctionReturnsComplex))
			{
				node->addNodeFlags(kAstReturnsComplex);
			}
		}


		bool allConst = false;
		bool allConstComplex = false;

		if (!sym ->hasSymbolFlag(kAstSymbolHasState))
		{
			allConst = true;
			allConstComplex = true;

			for (i = 0; i < count; i++)
			{
				auto n = node->getAt(i);

				allConst &= n->isImm() && !n->returnsComplex();
				allConstComplex &= allConst &= n->isImm() && n->returnsComplex();
			}
		}

		if (allConst && count <= 8) 
		{
			AstImm** args = reinterpret_cast<AstImm**>(node->getChildren());

			void* fn = sym->getFuncPtr();
			if (!sym->hasSymbolFlag(kAstSymbolRealFunctionReturnsComplex))
			{
				double result = 0.0;
#define ARG(n) args[n]->getValue()
				switch (count) 
				{
				case 0: result = ((Arg0Func)fn)(); break;
				case 1: result = ((Arg1Func)fn)(ARG(0)); break;
				case 2: result = ((Arg2Func)fn)(ARG(0), ARG(1)); break;
				case 3: result = ((Arg3Func)fn)(ARG(0), ARG(1), ARG(2)); break;
				case 4: result = ((Arg4Func)fn)(ARG(0), ARG(1), ARG(2), ARG(3)); break;
				case 5: result = ((Arg5Func)fn)(ARG(0), ARG(1), ARG(2), ARG(3), ARG(4)); break;
				case 6: result = ((Arg6Func)fn)(ARG(0), ARG(1), ARG(2), ARG(3), ARG(4), ARG(5)); break;
				case 7: result = ((Arg7Func)fn)(ARG(0), ARG(1), ARG(2), ARG(3), ARG(4), ARG(5), ARG(6)); break;
				case 8: result = ((Arg8Func)fn)(ARG(0), ARG(1), ARG(2), ARG(3), ARG(4), ARG(5), ARG(6), ARG(7)); break;
				}
#undef ARG

				AstNode* replacement = _ast->newNode<AstImm>(result);
				node->getParent()->replaceNode(node, replacement);
				onNode(replacement);
			}
			else
			{
				double argsDouble[9];
				for (i = 0; i < count; i++) {
					argsDouble[i] = static_cast<AstImm*>(node->getAt(i))->getValue();
				}
				std::complex<double> result(0.0, 0.0);
				result = ((mpFuncpDtoC)sym->getFuncPtr())(argsDouble);
				AstNode* replacement = _ast->newNode<AstImm>(result);
				node->getParent()->replaceNode(node, replacement);
				onNode(replacement);
			}
			_ast->deleteNode(node);
		}
		else if (allConstComplex && count < 9) 
		{
			std::complex<double> argsComplex[9];
			for (i = 0; i < count; i++) 
			{
				argsComplex[i] = static_cast<AstImm*>(node->getAt(i))->getValueComp();
			}

			AstImm* replacement = _ast->newNode<AstImm>(0);
			if (!sym->hasSymbolFlag(kAstSymbolComplexFunctionReturnsReal)) 
			{
				replacement->setValue(((mpFuncpCtoC)sym->getFuncPtr(true))(argsComplex));
			} 
			else 
			{
				replacement->setValue(((mpFuncpCtoD)sym->getFuncPtr(true))(argsComplex));
			}
			node->getParent()->replaceNode(node, replacement);
			onNode(replacement);
			
			_ast->deleteNode(node);
			
		}

		return kErrorOk;
	}



} // mathpresso namespace
