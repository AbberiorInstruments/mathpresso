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
#include "./mpoperation_p.h"

namespace mathpresso {

	// ============================================================================
	// [mpsl::AstOptimizer - Construction / Destruction]
	// ============================================================================

	AstOptimizer::AstOptimizer(AstBuilder* ast, ErrorReporter* errorReporter, const symbolMap* syms)
		: AstVisitor(ast),
		_errorReporter(errorReporter),
		_symbols(syms) {}
	AstOptimizer::~AstOptimizer() {}

	// ============================================================================
	// [mpsl::AstOptimizer - OnNode]
	// ============================================================================

	Error AstOptimizer::onBlock(AstBlock* node) {
		// Prevent removing nodes that are not stored in pure `AstBlock`. For example
		// function call inherits from `AstBlock`, but it needs each expression passed.
		bool alterable = node->getNodeType() == kAstNodeBlock;

		size_t curCount = node->getLength();
		size_t oldCount;
		bool isComplex = false;

		size_t i = 0;
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
				if (child->returnsComplex())
				{
					sym->setValue(static_cast<AstImm*>(child)->getValueCplx());
					node->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);
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

	Error AstOptimizer::onVar(AstVar* node) {
		AstSymbol* sym = node->getSymbol();
		bool b_complex = node->returnsComplex() || sym->hasSymbolFlag(kAstSymbolIsComplex);

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
		if (node->mpOp_)
		{
			node->mpOp_->optimize(this, node);
		}
		const OpInfo& op = OpInfo::get(node->getOp());

		MATHPRESSO_PROPAGATE(onNode(node->getChild()));
		AstNode* child = node->getChild();
		
		if (child->returnsComplex()) 
		{
			if (op.hasCtoC())
			{
				node->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);
			}
			else if (op.hasCtoD())
			{
				node->addNodeFlags(kAstTakesComplex);
			}
		}
		else 
		{
			if (op.hasDtoD()) 
			{
				//do nothing;
			} 
			else if (op.hasDtoC())
			{
				node->addNodeFlags(kAstReturnsComplex);
			}
			else if (op.hasCtoD()) 
			{
				node->addNodeFlags(kAstTakesComplex);
			}
			else if (op.hasCtoC()) 
			{
				node->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);
			}
		}

		if (child->isImm())
		{
			AstImm* p_imm = static_cast<AstImm*>(node->getChild());

			if (!p_imm ->returnsComplex())
			{// real parameter
				if (op.hasDtoD())
				{
					double value = p_imm->getValue();

					if (!op.isIntrinsic())
					{
						switch (node->getOp())
						{
						case kOpNeg: value = -value; break;
						case kOpNot: value = (value == 0); break;

						default:
							return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
								"Invalid unary operation '%s'.", op.name);
						}
					}
					else 
					{
						if (op.funcDtoD)
						{
							value = ((Arg1Func)op.funcDtoD)(value);
						}
						else 
						{
							return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
								"Invalid unary operation '%s'.", op.name);
						}
					}

					p_imm->setValue(value);
				}
				else if (op.hasDtoC()) 
				{
					double value = p_imm->getValue();
					std::complex<double> ret(0, 0);
					if (op.funcDtoD)
					{
						p_imm->setValue(((mpFuncpDtoC)op.funcDtoC)(&value));
					}
					else 
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid unary operation '%s'.", op.name);
					}
					node->addNodeFlags(kAstReturnsComplex);
				}
				else if (op.hasCtoD())
				{
					double operand = static_cast<AstImm*>(node->getChild())->getValue();
					std::complex<double> value(0, 0);
					if (op.funcCtoD)
					{
						p_imm->setValue(((mpFuncpCtoD)op.funcCtoD)(&std::complex<double>(operand, 0)));
					}
					else 
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"invalid unary operation '%s'.", op.name);
					}
				}
				else if (op.hasCtoC()) 
				{
					double operand = static_cast<AstImm*>(node->getChild())->getValue();
					std::complex<double> value(0, 0);
					if (op.funcCtoC) 
					{
						p_imm->setValue(((mpFuncpCtoC)op.funcCtoC)(&std::complex<double>(operand, 0)));
					}
					else 
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"invalid unary operation '%s'.", op.name);
					}
					node->addNodeFlags(kAstReturnsComplex);
				}
				else {
					return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
						"Invalid unary operation with complex parameters '%s'.", op.name);

				}

				node->unlinkChild();
				node->getParent()->replaceNode(node, p_imm);
				_ast->deleteNode(node);
			}
			else
			{ // Complex parameter
				std::complex<double> value = p_imm->getValueCplx();

				if (!op.isIntrinsic()) 
				{
					switch (node->getOp())
					{
					case kOpNeg: value = -value; break;
					case kOpNot: value = (value == std::complex<double>(0, 0)); break;

					default:
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid unary operation '%s'.", op.name);
					}
					p_imm->setValue(value);
				}
				else
				{
					if (op.hasCtoD())
					{
						p_imm->setValue(((mpFuncpCtoD)op.funcCtoD)(&value));
						p_imm->removeNodeFlags(kAstReturnsComplex);
					}
					else if (op.hasCtoC())
					{
						p_imm->setValue(((mpFuncpCtoC)op.funcCtoC)(&value));
					}
					else
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid unary operation with complex parameters '%s'.", op.name);
					}
				}
				node->unlinkChild();
				node->getParent()->replaceNode(node, p_imm);
				_ast->deleteNode(node);
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
			if (op.hasCtoD()) 
			{
				node->addNodeFlags(kAstTakesComplex);
			} 
			else if (op.hasCtoC()) 
			{
				node->addNodeFlags(kAstTakesComplex);
				node->addNodeFlags(kAstReturnsComplex);
			}
		}

		return kErrorOk;
	}

	Error AstOptimizer::onBinaryOp(AstBinaryOp* node)
	{
		const OpInfo& op = OpInfo::get(node->getOp());

		AstNode* left = node->getLeft();
		AstNode* right = node->getRight();

		if (node->mpOp_ != nullptr)
		{
			return node->mpOp_->optimize(this, node);
		}

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

			/*if (!OpInfo::get(branchCondition->getOp()).isCondition()) {
				return _errorReporter->onError(kErrorInvalidArgument, node->getPosition(),
					"Not a condition.");
			}*/
			// remove branchCondition from the AST
			node->setLeft(nullptr);
			branchCondition->_parent = nullptr;

			// remove the right path from the AST.
			lastColon->setRight(nullptr);
			branchRight->_parent = nullptr;


			// Distinguish between a complex and a non-complex case:
			// i.e.: cond1 ? cond2 ? a : b : c
			if (node->getRight() != lastColon) {
				// remove left branch from the AST.
				branchLeft = node->getRight();
				node->setRight(nullptr);
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
			AstTernaryOp* newNode = node->getAst()->newNode<AstTernaryOp>(kOpQMark);
			newNode->setCondition(branchCondition);
			newNode->setLeft(branchLeft);
			newNode->setRight(branchRight);

			AstBinaryOp* oldNode = node;

			// add the new node to the AST.
			node->getParent()->replaceNode(node, newNode);

			// clean up:
			lastColon->setLeft(nullptr);
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
			// If both nodes are values it's easy, just fold them into a single one.
			if (lIsImm && rIsImm)
			{
				AstImm* lNode = static_cast<AstImm*>(left);
				AstImm* rNode = static_cast<AstImm*>(right);

				double lVal = lNode->getValue();
				double rVal = rNode->getValue();
				if (op.hasDtoD())
				{
					if (!op.isIntrinsic())
					{
						if (op.funcDtoD) 
						{
							lNode->setValue(((Arg2Func)op.funcDtoD)(lVal, rVal));
						}
						else
						{
							double result = 0.0;
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
							lNode->setValue(result);
						}
					}
					else
					{
						if (op.funcDtoD)
						{
							lNode->setValue(((Arg2Func)op.funcDtoD)(lVal, rVal));
						}
						else
						{
							return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
								"Invalid binary operation '%s'.", op.name);
						}
					}
				}
				else if (op.hasDtoC())
				{
					double args[] = { lVal, rVal };
					if (op.funcDtoD)
					{
						lNode->setValue(((mpFuncpDtoC)op.funcDtoC)(args));
					}
					else
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid binary operation '%s'.", op.name);
					}
				}
				else if (op.hasCtoD())
				{
					std::complex<double> args[] = { {lVal,0}, {rVal,0} };
					if (op.funcDtoD)
					{
						lNode->setValue(((mpFuncpCtoD)op.funcCtoD)(args));
					}
					else
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid binary operation '%s'.", op.name);
					}
				}
				else if (op.hasCtoC())
				{
					std::complex<double> args[] = { { lVal,0 },{ rVal,0 } };
					if (op.funcDtoD)
					{
						lNode->setValue(((mpFuncpCtoC)op.funcCtoD)(args));
					}
					else
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid binary operation '%s'.", op.name);
					}
				}
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
			if (!op.hasCtoD() && op.hasCtoC()) 
			{
				node->addNodeFlags(kAstReturnsComplex);
			}
			
			// if we have to calculate in complex, and one of the operands is an immediate, it should be converted to complex.
			if (left->isImm()&& !left->returnsComplex()) 
			{
				left->addNodeFlags(kAstReturnsComplex);
			} 
			else if (right->isImm() && !right->returnsComplex())
			{
				right->addNodeFlags(kAstReturnsComplex);
			}

			// complex immediates.
			if (left->isImm() && right->isImm() && left ->returnsComplex() && right ->returnsComplex()) // both should be complex, but who cares ;)
			{
				AstImm* lNode = static_cast<AstImm*>(left);
				AstImm* rNode = static_cast<AstImm*>(right);

				if (!op.isIntrinsic()) {
					std::complex<double> lVal = lNode->getValueCplx();
					std::complex<double> rVal = rNode->getValueCplx();
					if (op.funcDtoD)
					{
						std::complex<double> args[] = { lVal, rVal};
						rNode->setValue(((mpFuncpCtoC)op.funcCtoC)(args));
					}
					else
					{
						std::complex<double> result;
						switch (node->getOp())
						{
						case kOpAdd: result = lVal + rVal; break;
						case kOpSub: result = lVal - rVal; break;
						case kOpMul: result = lVal * rVal; break;
						case kOpDiv: result = lVal / rVal; break;
						case kOpEq: result = lVal == rVal; break;
						case kOpNe: result = lVal != rVal; break;

						default:
							return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
								"Invalid complex binary operation '%s'.", op.name);
						}
						rNode->setValue(result);
					}
				}
				else if (op.hasCtoD()) 
				{
					std::complex<double> args[] = { lNode->getValueCplx(), rNode->getValueCplx() };
					if (op.funcCtoD) 
					{
						rNode->setValue(((mpFuncpCtoD)op.funcCtoD)(args));
					}
					else 
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid complex binary operation '%s'.", op.name);
					}
				} 
				else  if (op.hasCtoC())
				{
					std::complex<double> args[] = { lNode->getValueCplx(), rNode->getValueCplx() };
					if (op.funcCtoC) 
					{
						rNode->setValue(((mpFuncpCtoC)op.funcCtoC)(args));
					}
					else 
					{
						return _errorReporter->onError(kErrorInvalidState, node->getPosition(),
							"Invalid complex binary operation '%s'.", op.name);
					}
				}

				node->unlinkRight();
				node->getParent()->replaceNode(node, rNode);

				_ast->deleteNode(node);
				onNode(rNode);
				return kErrorOk;
			}
			else if (lIsImm)
			{
				AstImm* lNode = static_cast<AstImm*>(left);
				std::complex<double> val = lNode->getValueCplx();

				if ((val == std::complex<double>(0.0, 0.0) && (op.flags & kOpFlagNopIfLZero)) ||
					(val == std::complex<double>(1.0, 0.0) && (op.flags & kOpFlagNopIfLOne)))
				{
					node->unlinkRight();
					node->getParent()->replaceNode(node, right);

					_ast->deleteNode(node);
				}
			}
			else if (rIsImm) 
			{
				AstImm* rNode = static_cast<AstImm*>(right);
				std::complex<double> val = rNode->getValueCplx();

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
					if ((val == std::complex<double>(0.0, 0.0) && (op.flags & kOpFlagNopIfRZero)) ||
						(val == std::complex<double>(1.0, 0.0) && (op.flags & kOpFlagNopIfROne)))
					{
						node->unlinkLeft();
						node->getParent()->replaceNode(node, left);

						_ast->deleteNode(node);
					}
				}
			}

		}

		return kErrorOk;
	}

	Error AstOptimizer::onTernaryOp(AstTernaryOp* node) {

		if (node->mpOp_ != nullptr)
		{
			return node->mpOp_->optimize(this, node);
		}

		MATHPRESSO_PROPAGATE(onNode(node->getCondition()));
		AstNode* branchCond = node->getCondition();
		if (branchCond->isImm()) {
			bool conditionIsTrue = static_cast<AstImm*>(branchCond)->getValueCplx() != std::complex<double>({ 0, 0 });
			
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

	Error AstOptimizer::onCall(AstCall* node) 
	{
		return node->getSymbol()->getOp()->optimize(this, node);
	}



} // mathpresso namespace
