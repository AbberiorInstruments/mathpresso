// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include <mathpresso/mpparser_p.h>
#include <mathpresso/mpoperation.h>

namespace mathpresso {

	// ============================================================================
	// [mathpresso::Parser - Syntax Error]
	// ============================================================================

#define MATHPRESSO_PARSER_ERROR(_Token_, ...) \
  return _errorReporter->onError( \
    ErrorCode::kErrorInvalidSyntax, static_cast<uint32_t>((size_t)(_Token_.position)), __VA_ARGS__)

#define MATHPRESSO_PARSER_WARNING(_Token_, ...) \
  _errorReporter->onWarning( \
    static_cast<uint32_t>((size_t)(_Token_.position)), __VA_ARGS__)

// ============================================================================
// [mathpresso::AstNestedScope]
// ============================================================================

//! \internal
//!
//! Nested scope used only by the parser and always allocated statically.
	struct AstNestedScope : public AstScope {
		MATHPRESSO_NO_COPY(AstNestedScope)

			// --------------------------------------------------------------------------
			// [Construction / Destruction]
			// --------------------------------------------------------------------------

			MATHPRESSO_INLINE AstNestedScope(Parser* parser)
			: AstScope(parser->_ast, parser->_currentScope, AstScopeType::kAstScopeNested),
			_parser(parser)
		{
			_parser->_currentScope = this;
		}

		MATHPRESSO_INLINE ~AstNestedScope()
		{
			AstScope* parent = getParent();
			MATHPRESSO_ASSERT(parent != nullptr);

			_parser->_currentScope = parent;
			parent->_symbols.mergeToInvisibleSlot(this->_symbols);
		}

		// --------------------------------------------------------------------------
		// [Members]
		// --------------------------------------------------------------------------

		Parser* _parser;
	};

	// ============================================================================
	// [mathpresso::Parser - Parse]
	// ============================================================================

	Error Parser::parseProgram(AstProgram* block)
	{
		for (;;)
		{
			Token token;
			uint32_t uToken = _tokenizer.peek(&token);

			// Parse the end of the input.
			if (uToken == TokenType::kTokenEnd)
				break;
			//MATHPRESSO_PROPAGATE(parseStatement(block, ParserFlags::kEnableVarDecls | ParserFlags::kEnableNestedBlock));
			auto ret = parseStatement(block, ParserFlags::kEnableVarDecls | ParserFlags::kEnableNestedBlock);
			if (ret == ErrorCode::kErrorOk)
			{
				reparseTernary(block);
			}
			else
			{
				return ret;
			}
		}

		if (block->getLength() == 0)
			return MATHPRESSO_TRACE_ERROR(ErrorCode::kErrorNoExpression);

		return ErrorCode::kErrorOk;
	}

	// Parse <statement>; or { [<statement>; ...] }
	Error Parser::parseStatement(AstBlock* block, uint32_t flags)
	{
		Token token;
		uint32_t uToken = _tokenizer.peek(&token);

		// Parse the ';' token.
		if (uToken == TokenType::kTokenSemicolon)
		{
			_tokenizer.consume();
			return ErrorCode::kErrorOk;
		}


		// Parse a nested block.
		if (uToken == TokenType::kTokenLCurl)
		{
			AstBlock* nested;

			if (!(flags & ParserFlags::kEnableNestedBlock))
				MATHPRESSO_PARSER_ERROR(token, "Cannot declare a new block-scope here.");

			MATHPRESSO_PROPAGATE(block->willAdd());
			MATHPRESSO_NULLCHECK(nested = _ast->newNode<AstBlock>());
			block->appendNode(nested);

			AstNestedScope tmpScope(this);
			return parseBlockOrStatement(nested);
		}

		// Parse a variable declaration.
		if (uToken == TokenType::kTokenVar)
		{
			if (!(flags & ParserFlags::kEnableVarDecls))
				MATHPRESSO_PARSER_ERROR(token, "Cannot declare a new variable here.");
			return parseVariableDecl(block);
		}

		// Parse an expression.
		AstNode* expression;

		MATHPRESSO_PROPAGATE(block->willAdd());
		MATHPRESSO_PROPAGATE(parseExpression(&expression, false));
		block->appendNode(expression);

		uToken = _tokenizer.peek(&token);
		if (uToken == TokenType::kTokenSemicolon)
		{
			_tokenizer.consume();
			return ErrorCode::kErrorOk;
		}

		if (uToken == TokenType::kTokenEnd)
			return ErrorCode::kErrorOk;

		MATHPRESSO_PARSER_ERROR(token, "Expected a ';' after an expression.");
	}

	// Parse <block|statement>;.
	Error Parser::parseBlockOrStatement(AstBlock* block)
	{
		Token token;
		uint32_t uToken = _tokenizer.next(&token);

		// Parse the <block>, consume '{' token.
		block->setPosition(token.getPosAsUInt());
		if (uToken == TokenType::kTokenLCurl)
		{
			for (;;)
			{
				uToken = _tokenizer.peek(&token);

				// Parse the end of the block '}'.
				if (uToken == TokenType::kTokenRCurl)
				{
					// Consume the '}' token.
					_tokenizer.consume();
					return ErrorCode::kErrorOk;
				}
				else
				{
					MATHPRESSO_PROPAGATE(parseStatement(block, ParserFlags::kEnableVarDecls | ParserFlags::kEnableNestedBlock));
				}
			}
		}
		else
		{
			return parseStatement(block, ParserFlags::kNoFlags);
			
		}
	}

	// Parse "var <name> = <expression>[, <name> = <expression>, ...];".
	Error Parser::parseVariableDecl(AstBlock* block)
	{
		Token token;
		uint32_t uToken = _tokenizer.next(&token);
		std::string str;

		bool isFirst = true;
		uint32_t position = token.getPosAsUInt();

		// Parse the 'var' keyword.
		if (uToken != TokenType::kTokenVar)
			MATHPRESSO_PARSER_ERROR(token, "Expected 'var' keyword.");

		AstScope* scope = _currentScope;
		for (;;)
		{
			// Parse the variable name.
			if (_tokenizer.next(&token) != TokenType::kTokenSymbol)
				MATHPRESSO_PARSER_ERROR(token, isFirst
					? "Expected a variable name after 'var' keyword."
					: "Expected a variable name after colon ','.");

			MATHPRESSO_PROPAGATE(block->willAdd());

			if (!isFirst)
				position = token.getPosAsUInt();

			// Resolve the variable name.
			AstSymbol* vSym;
			AstScope* vScope;

			str = std::string(_tokenizer._start + token.position, token.length);
			if ((vSym = scope->resolveSymbol(str, token.hVal, &vScope)) != nullptr)
			{
				if (vSym->getSymbolType() != AstSymbolType::kAstSymbolVariable || scope == vScope)
					MATHPRESSO_PARSER_ERROR(token, "Attempt to redefine '%s'.", vSym->getName());

				if (vSym->hasNode())
				{
					uint32_t line, column;
					_errorReporter->getLineAndColumn(vSym->getNode()->getPosition(), line, column);
					MATHPRESSO_PARSER_WARNING(token, "Variable '%s' shadows a variable declared at [%d:%d].", vSym->getName(), line, column);
				}
				else
				{
					MATHPRESSO_PARSER_WARNING(token, "Variable '%s' shadows a variable of the same name.", vSym->getName());
				}
			}

			vSym = _ast->newSymbol(str, token.hVal, AstSymbolType::kAstSymbolVariable, scope->getScopeType());
			MATHPRESSO_NULLCHECK(vSym);
			scope->putSymbol(vSym);

			AstVarDecl* decl = _ast->newNode<AstVarDecl>();
			MATHPRESSO_NULLCHECK_(decl, { _ast->deleteSymbol(vSym); });
			decl->_mpOp = _ops->find("=", 2).get();

			decl->setPosition(position);
			decl->setSymbol(vSym);

			// Assign a slot and fill to safe defaults.
			vSym->setVarOffset(0);
			vSym->setVarSlotId(_ast->newSlotId());
			vSym->setNode(decl);

			// Parse possible assignment '='.
			uToken = _tokenizer.next(&token);
			bool isAssigned = (uToken == TokenType::kTokenOperator && std::string(_tokenizer._start + token.position, token.length) == "=");

			if (isAssigned)
			{
				AstNode* expression;
				MATHPRESSO_PROPAGATE_(parseExpression(&expression, false), { _ast->deleteNode(decl); });

				decl->setChild(expression);
				vSym->incWriteCount();

				uToken = _tokenizer.next(&token);
			}

			// Make the symbol declared so it can be referenced after now.
			vSym->setDeclared();

			// Parse the ',' or ';' tokens.
			if (uToken == TokenType::kTokenComma || uToken == TokenType::kTokenSemicolon || uToken == TokenType::kTokenEnd)
			{
				block->appendNode(decl);

				// Token ';' terminates the declaration.
				if (uToken != TokenType::kTokenComma)
					break;
			}
			else
			{
				_ast->deleteSymbol(vSym);
				MATHPRESSO_PARSER_ERROR(token, "Unexpected token %d.", token.token);
			}

			isFirst = false;
		}

		return ErrorCode::kErrorOk;
	}

	Error Parser::parseExpression(AstNode** pNode, bool isNested)
	{
		AstScope* scope = _currentScope;

		// It's important that the given expression is parsed in a way that it can be
		// correctly evaluated. The `parseExpression()` function can handle expressions
		// that contain lastUnaryNode and binary operators combined with terminals (variables,
		// constants or function calls).
		//
		// The most expression parsers usually use stack to handle operator precedence,
		// but MPSL uses AstNode parent->child hierarchy instead. When operator with a
		// higher precedence is found it traverses down in the hierarchy and when
		// operator with the same/lower precedence is found the hierarchy is traversed
		// back.
		//
		//                               AST Examples
		//
		// +-----------------+-----------------+-----------------+-----------------+
		// |                 |                 |                 |                 |
		// |   "a + b + c"   |   "a * b + c"   |   "a + b * c"   |   "a * b * c"   |
		// |                 |                 |                 |                 |
		// |       (+)       |       (+)       |       (+)       |       (*)       |
		// |      /   \      |      /   \      |      /   \      |      /   \      |
		// |   (+)     (c)   |   (*)     (c)   |   (a)     (*)   |   (*)     (c)   |
		// |   / \           |   / \           |           / \   |   / \           |
		// | (a) (b)         | (a) (b)         |         (b) (c) | (a) (b)         |
		// |                 |                 |                 |                 |
		// +-----------------+-----------------+-----------------+-----------------+
		// |                 |                 |                 |                 |
		// | "a + b + c + d" | "a + b * c + d" | "a * b + c * d" | "a * b * c + d" |
		// |                 |                 |                 |                 |
		// |       (+)       |       (+)       |       (+)       |       (+)       |
		// |      /   \      |      /   \      |      /   \      |      /   \      |
		// |   (+)     (d)   |    (+)   (d)    |   (*)     (*)   |   (*)     (d)   |
		// |    | \          |   /   \         |   / \     / \   |    | \          |
		// |   (+) (c)       |(a)     (*)      | (a) (b) (c) (d) |   (*) (c)       |
		// |   / \           |        / \      |                 |   / \           |
		// | (a) (b)         |      (b) (c)    |                 | (a) (b)         |
		// |                 |                 |                 |                 |
		// +-----------------+-----------------+-----------------+-----------------+

		Token token;
		
		// Current binary operator node. Initial nullptr value means that the parsing
		// just started and there is no binary operator yet. Once the first binary
		// operator has been parsed `currentBinaryNode` will be set accordingly.
		AstBinaryOp* currentBinaryNode = nullptr;

		// Currently parsed node.
		AstNode* currentNode = nullptr;

		for (;;)
		{
			// Last unary node. It's an optimization to prevent recursion in case that
			// we found two or more unary expressions after each other. For example the
			// expression "-!-1" contains only unary operators that will be parsed by
			// a single `parseExpression()` call.
			AstUnary* lastUnaryNode = nullptr;
			bool b_complex = false;

		_Repeat1:
			switch (_tokenizer.next(&token))
			{
				// Parse a variable, a constant, or a function-call. This can be repeated
				// one or several times based on the expression type. For unary nodes it's
				// repeated immediately, for binary nodes it's repeated after the binary
				// node is created.

				// Parse a symbol (variable or function name).
			case TokenType::kTokenSymbol:
			{
				std::string str(_tokenizer._start + token.position, token.length);

				AstScope* symScope;
				AstSymbol* sym = scope->resolveSymbol(str, token.hVal, &symScope);

				if (sym == nullptr)
					MATHPRESSO_PARSER_ERROR(token, "Unresolved symbol %.*s.", static_cast<int>(str.length()), str.c_str());

				uint32_t symType = sym->getSymbolType();
				AstNode* newNode;

				if (symType == AstSymbolType::kAstSymbolVariable)
				{
					if (!sym->isDeclared())
						MATHPRESSO_PARSER_ERROR(token, "Can't use variable '%s' that is being declared.", sym->getName());

					// Put symbol to shadow scope if it's global. This is done lazily and
					// only once per symbol when it's referenced.
					if (symScope->isGlobal())
					{
						sym = _ast->shadowSymbol(sym);
						MATHPRESSO_NULLCHECK(sym);

						sym->setVarSlotId(_ast->newSlotId());
						symScope = _ast->getRootScope();
						symScope->putSymbol(sym);
					}

					newNode = _ast->newNode<AstVar>();
					MATHPRESSO_NULLCHECK(newNode);
					static_cast<AstVar*>(newNode)->setSymbol(sym);

					if (sym->hasSymbolFlag(AstSymbolFlags::kAstSymbolIsComplex))
						newNode->addNodeFlags(AstNodeFlags::kAstTakesComplex | AstNodeFlags::kAstReturnsComplex);

					newNode->setPosition(token.getPosAsUInt());
					sym->incUsedCount();
				}
				else
				{
					// Will be parsed by `parseCall()` again.
					_tokenizer.set(&token);
					MATHPRESSO_PROPAGATE(parseCall(&newNode));
				}

				if (lastUnaryNode == nullptr)
					currentNode = newNode;
				else
					lastUnaryNode->setChild(newNode);
				break;
			}

			// Parse a number.
			case TokenType::kTokenComplex:
				b_complex = true;
			case TokenType::kTokenNumber:
			{
				AstImm* newNode = _ast->newNode<AstImm>();
				MATHPRESSO_NULLCHECK(newNode);

				newNode->setPosition(token.getPosAsUInt());
				if (!b_complex)
				{
					newNode->setValue(token.value);
				}
				else
				{
					newNode->setValue({ 0, token.value });
				}

				if (lastUnaryNode == nullptr)
					currentNode = newNode;
				else
					lastUnaryNode->setChild(newNode);
				break;
			}

			// Parse expression terminators - ',', ';' or ')'.
			case TokenType::kTokenComma:
			case TokenType::kTokenSemicolon:
			case TokenType::kTokenRParen:
			{
				MATHPRESSO_PARSER_ERROR(token, "Expected an expression.");
			}

			// Parse a nested expression.
			case TokenType::kTokenLParen:
			{
				uint32_t position = token.getPosAsUInt();

				AstNode* newNode;
				MATHPRESSO_PROPAGATE(parseExpression(&newNode, true));

				if (_tokenizer.next(&token) != TokenType::kTokenRParen)
					MATHPRESSO_PARSER_ERROR(token, "Expected a ')' token.");

				if (lastUnaryNode == nullptr)
					currentNode = newNode;
				else
					lastUnaryNode->setChild(newNode);

				break;
			}

			// Parse a right-to-left associative unary operator ('+', '-', "!").
			case TokenType::kTokenOperator:
			{
				std::string name(_tokenizer._start + token.position, token.length);
				auto op = _ops->find(name, 1);
				if (!op)
					MATHPRESSO_PARSER_ERROR(token, "Invalid unary operator.");

				
				// Parse the unary operator.
				AstUnaryOp* opNode = _ast->newNode<AstUnaryOp>();
				MATHPRESSO_NULLCHECK(opNode);
				opNode->setPosition(token.getPosAsUInt());

				opNode->_mpOp = op.get();
				
				if (lastUnaryNode == nullptr)
					currentNode = opNode;
				else
					lastUnaryNode->setChild(opNode);

				isNested = true;
				lastUnaryNode = opNode;

				goto _Repeat1;
			}

			case TokenType::kTokenEnd:
			{
				MATHPRESSO_PARSER_ERROR(token, "Unexpected end of the program.");
			}

			default:
			{
				MATHPRESSO_PARSER_ERROR(token, "Unexpected token %d.", token.token);
			}
			}

			// _Repeat2:
			switch (_tokenizer.next(&token))
			{
				// Parse the expression terminators - ',', ';', ')' or EOI.
			case TokenType::kTokenComma:
			case TokenType::kTokenSemicolon:
			case TokenType::kTokenRParen:
			case TokenType::kTokenEnd:
			{
				_tokenizer.set(&token);

				if (currentBinaryNode != nullptr)
				{
					currentBinaryNode->setRight(currentNode);
					// Iterate to the top-most node.
					while (currentBinaryNode->hasParent())
						currentBinaryNode = static_cast<AstBinaryOp*>(currentBinaryNode->getParent());
					currentNode = currentBinaryNode;
				}

				*pNode = currentNode;
				return ErrorCode::kErrorOk;
			}

			// Parse Binary Operators
			case TokenType::kTokenOperator:
			{
				std::string name(_tokenizer._start + token.position, token.length);
				auto op = _ops->find(name, 2);
				if (!op)
					MATHPRESSO_PARSER_ERROR(token, "Invalid Operator.");


				if (name == "=")
				{
					// Check whether the assignment is valid.
					if (currentNode->getNodeType() != AstNodeType::kAstNodeVar)
						MATHPRESSO_PARSER_ERROR(token, "Can't assign to a non-variable.");

					AstSymbol* sym = static_cast<AstVar*>(currentNode)->getSymbol();
					if (sym->hasSymbolFlag(AstSymbolFlags::kAstSymbolIsReadOnly))
						MATHPRESSO_PARSER_ERROR(token, "Can't assign to a read-only variable '%s'.", sym->getName());

					if (isNested)
						MATHPRESSO_PARSER_ERROR(token, "Invalid assignment inside an expression.");

					sym->incWriteCount();
				}

				AstBinaryOp* newNode = _ast->newNode<AstBinaryOp>();
				MATHPRESSO_NULLCHECK(newNode);

				newNode->_mpOp = op.get();

				newNode->setPosition(token.getPosAsUInt());

				if (currentBinaryNode == nullptr)
				{
					// currentBinaryNode <------+
					//              |
					// +------------+------------+ First operand - currentBinaryNode becomes the newly
					// |        (newNode)        | created newNode; currentNode is assigned to the
					// |        /       \        | left side of newNode and will be referred
					// | (currentNode)  (NULL)   | as (...) by the next operation.
					// +-------------------------+
					newNode->setLeft(currentNode);
					currentBinaryNode = newNode;
					break;
				}

				uint32_t currentBinaryPrec = currentBinaryNode->_mpOp->precedence();
				uint32_t newBinaryPrec = newNode->_mpOp->precedence();

				if (currentBinaryPrec > newBinaryPrec)
				{
					// currentBinaryNode <-+
					//                     |
					// +-------------------+-----+ The current operator (newBinaryPrec) has a
					// |(currentBinaryNode)|     | higher precedence than the previous one
					// |    /       \      |     | (currentBinaryPrec), so the newNode will be assigned
					// | (...)      (newNode)    | to the right side of currentBinaryNode and it will
					// |            /       \    | function as a stack-like structure. We
					// |    (currentNode)  (NULL)| have to advance back at some point.
					// +-------------------------+

					currentBinaryNode->setRight(newNode);
					newNode->setLeft(currentNode);
					currentBinaryNode = newNode;
					break;
				}
				else
				{
					currentBinaryNode->setRight(currentNode);

					// Advance to the top-most binaryNode that has less or equal precedence
					// than newBinaryPrec.
					while (currentBinaryNode->hasParent())
					{
						// Terminate conditions:
						//   1. currentBinaryNode has higher precedence than newNode.
						//   2. currentBinaryNode has equal precedence and right-to-left associativity.
						if (currentBinaryPrec > newBinaryPrec || (currentBinaryPrec == newBinaryPrec && currentBinaryNode->_mpOp->isRightToLeft()))
							break;
						currentBinaryNode = static_cast<AstBinaryOp*>(currentBinaryNode->getParent());
						currentBinaryPrec = currentBinaryNode->_mpOp->precedence();
					}

					// currentBinaryNode <+
					//                    |
					// +------------------+------+
					// |           (newNode)     | Simple case - currentBinaryNode becomes the left
					// |           /       \     | node in the created newNode and newNode
					// |(currentBinaryNode)(NULL)| becomes currentBinaryNode for the next operator.
					// |    /       \            |
					// | (...)    (currentNode)  | currentBinaryNode will become a top-level node.
					// +-------------------------+

					if (!currentBinaryNode->hasParent() && 
						!(currentBinaryPrec > newBinaryPrec || (currentBinaryPrec == newBinaryPrec && currentBinaryNode->_mpOp->isRightToLeft())))
					{
						newNode->setLeft(currentBinaryNode);
					}
					// currentBinaryNode <-+
					//                     |
					// +-------------------+-----+
					// |(currentBinaryNode)|     |
					// |    /       \      |     | Complex case - inject node in place
					// | (...)      (newNode)    | of currentBinaryNode.right (because of higher
					// |            /       \    | precedence or RTL associativity).
					// |  (currentNode)   (NULL) |
					// +-------------------------+
					else
					{
						AstNode* pNode = currentBinaryNode->unlinkRight();
						currentBinaryNode->setRight(newNode);
						newNode->setLeft(pNode);
					}

					isNested = true;
					currentBinaryNode = newNode;

					break;
				}
			}



			default:
			{
				MATHPRESSO_PARSER_ERROR(token, "Unexpected token %d.", token.token);
			}
			}
		}
	}

	// Parse "function([arg1 [, arg2, ...] ])".
	Error Parser::parseCall(AstNode** pNodeOut)
	{
		Token token;
		uint32_t uToken;

		uToken = _tokenizer.next(&token);
		MATHPRESSO_ASSERT(uToken == TokenType::kTokenSymbol);
		uint32_t position = token.getPosAsUInt();

		std::string fnName(_tokenizer.getTokenName());

		uToken = _tokenizer.next(&token);
		if (uToken != TokenType::kTokenLParen)
			MATHPRESSO_PARSER_ERROR(token, "Expected a '(' token after a function name."); 

		AstCall* callNode = _ast->newNode<AstCall>(); 
		MATHPRESSO_NULLCHECK(callNode);

		callNode->setPosition(position);

		uToken = _tokenizer.peek(&token);
		if (uToken != TokenType::kTokenRParen)
		{ 
			for (;;)
			{
				// Parse the argument expression.
				AstNode* expression;
				Error err;

				if ((err = callNode->willAdd()) != ErrorCode::kErrorOk || (err = parseExpression(&expression, true)) != ErrorCode::kErrorOk)
				{
					_ast->deleteNode(callNode);
					return err;
				}

				callNode->appendNode(expression);

				// Parse ')' or ',' tokens.
				uToken = _tokenizer.peek(&token);
				if (uToken == TokenType::kTokenRParen)
					break;

				if (uToken != TokenType::kTokenComma)
				{
					_ast->deleteNode(callNode);
					MATHPRESSO_PARSER_ERROR(token, "Expected either ',' or ')' token.");
				}

				_tokenizer.consume();
			}
		}

		_tokenizer.consume();

		if (_ops->find(fnName, callNode->getLength()))
		{
			callNode->_mpOp = _ops->find(fnName, callNode->getLength()).get();
		}
		else
		{
			MATHPRESSO_PARSER_ERROR(token, "Function '%s' requires a MpOperation-Object with %d arguments.", fnName, callNode->getLength());
		}

		*pNodeOut = callNode;
		return ErrorCode::kErrorOk;
	}

	Error Parser::reparseTernary(AstNode * node)
	{

		for (size_t i = 0; i < node->getLength(); i++)
		{
			unsigned int ret;
			if (node->getNodeType() == AstNodeType::kAstNodeBinaryOp && _ops->name(node->_mpOp) == "?")
			{
				
				AstBinaryOp* lastColon = static_cast<AstBinaryOp*>(node);
				// go to the last Colon after question-marks.
				while (lastColon->_mpOp && _ops->name(lastColon->_mpOp) == "?")
				{
					lastColon = static_cast<AstBinaryOp*>(lastColon->getRight());
				}

				while (lastColon->getRight() && lastColon->getRight()->_mpOp && _ops->name(lastColon->getRight()->_mpOp) == ":") // check whether colon or Qmark
				{
					lastColon = static_cast<AstBinaryOp*>(lastColon->getRight());
				}

				if (_ops->name(lastColon->_mpOp) != ":")
				{
					return _errorReporter->onError(ErrorCode::kErrorInvalidSyntax, node->getPosition(),
														"Invalid ternary operation. Expected a ':'.");
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
				AstTernaryOp* ternaryNode = node->getAst()->newNode<AstTernaryOp>();
				ternaryNode->setCondition(branchCondition);
				ternaryNode->setLeft(branchLeft);
				ternaryNode->setRight(branchRight);
				ternaryNode->_mpOp = _ops->find("?", 2).get();

				// add the new node to the AST.
				node->getParent()->replaceNode(node, ternaryNode);

				// clean up:
				lastColon->setLeft(nullptr);
				_ast->deleteNode(lastColon);
				_ast->deleteNode(node);
				return reparseTernary(ternaryNode);
			}
			else
			{
				if (node->getAt(i))
				{
					ret = reparseTernary(node->getAt(i));
				}
				else
				{
					ret = ErrorCode::kErrorOk;
				}
			}

			if (ret != ErrorCode::kErrorOk)
			{
				return ret;
			}

		}

		return ErrorCode::kErrorOk;
	}

} // mathpresso namespace
