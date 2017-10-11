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

namespace mathpresso
{

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
	// [mathpresso::Parser - Parse]
	// ============================================================================

	Error Parser::parseProgram(std::shared_ptr<AstProgram> block)
	{
		for (;;)
		{
			Token token;
			uint32_t uToken = _tokenizer.peek(&token);

			// Parse the end of the input.
			if (uToken == TokenType::kTokenEnd)
				break;
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
	Error Parser::parseStatement(std::shared_ptr<AstBlock> block, uint32_t flags)
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
			std::shared_ptr<AstBlock> nested;

			if (!(flags & ParserFlags::kEnableNestedBlock))
				MATHPRESSO_PARSER_ERROR(token, "Cannot declare a new block-scope here.");

			MATHPRESSO_NULLCHECK(nested = std::make_shared<AstBlock>());
			block->appendNode(nested);

			auto nestedContext(std::make_shared<Context>());

			nestedContext->markShadow();
			nestedContext->setParent(_shadowContext);

			_shadowContext = nestedContext;

			auto ret = parseBlockOrStatement(nested);

			_shadowContext = nestedContext->getParent();

			return ret;
		}

		// Parse a variable declaration.
		if (uToken == TokenType::kTokenVar)
		{
			if (!(flags & ParserFlags::kEnableVarDecls))
				MATHPRESSO_PARSER_ERROR(token, "Cannot declare a new variable here.");
			return parseVariableDecl(block);
		}

		// Parse an expression.
		std::shared_ptr<AstNode> expression;

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
	Error Parser::parseBlockOrStatement(std::shared_ptr<AstBlock> block)
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
	Error Parser::parseVariableDecl(std::shared_ptr<AstBlock> block)
	{
		Token token;
		uint32_t uToken = _tokenizer.next(&token);
		std::string symbolName;

		bool isFirst = true;
		uint32_t position = token.getPosAsUInt();

		// Parse the 'var' keyword.
		if (uToken != TokenType::kTokenVar)
			MATHPRESSO_PARSER_ERROR(token, "Expected 'var' keyword.");

		for (;;)
		{
			// Parse the variable name.
			if (_tokenizer.next(&token) != TokenType::kTokenSymbol)
				MATHPRESSO_PARSER_ERROR(token, isFirst
										? "Expected a variable name after 'var' keyword."
										: "Expected a variable name after colon ','.");

			if (!isFirst)
				position = token.getPosAsUInt();

			// Resolve the variable name.
			std::shared_ptr<AstSymbol> vSym;
			std::shared_ptr<Context> vContext;

			symbolName = _tokenizer.getTokenName(token);
			if ((vSym = resolver::resolveVariable(_shadowContext, symbolName, &vContext)))
			{
				if (vSym->getSymbolType() != AstSymbolType::kAstSymbolVariable || _shadowContext == vContext)
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

			vSym = std::make_shared<AstSymbol>(symbolName, AstSymbolType::kAstSymbolVariable, _shadowContext->isGlobal());
			_shadowContext->_symbols->add(symbolName, vSym);
			MATHPRESSO_NULLCHECK(vSym);

			std::shared_ptr<AstVarDecl> decl = std::make_shared<AstVarDecl>();
			MATHPRESSO_NULLCHECK(decl);
			decl->_mpOp = resolver::resolveFunction(_shadowContext, "=", 2);
			decl->_opName = "=";

			decl->setPosition(position);
			decl->setSymbol(vSym);

			// Assign a slot and fill to safe defaults.
			vSym->setVarOffset(0);
			vSym->setVarSlotId(_ast->newSlotId());
			vSym->setNode(decl);

			// Parse possible assignment '='.
			uToken = _tokenizer.next(&token);
			bool isAssigned = (uToken == TokenType::kTokenOperator && _tokenizer.getTokenName(token) == "=");

			if (isAssigned)
			{
				std::shared_ptr<AstNode> expression;
				MATHPRESSO_PROPAGATE(parseExpression(&expression, false));

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
				MATHPRESSO_PARSER_ERROR(token, "Unexpected token %d.", token.token);
			}

			isFirst = false;
		}

		return ErrorCode::kErrorOk;
	}

	Error Parser::parseExpression(std::shared_ptr<AstNode>* pNodeOut, bool isNested)
	{

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
		std::shared_ptr<AstBinaryOp> currentBinaryNode = nullptr;

		// necessary to not loose 
		std::shared_ptr<AstBinaryOp> rootBinaryNode = nullptr;

		// Currently parsed node.
		std::shared_ptr<AstNode> currentNode = nullptr;

		for (;;)
		{
			// Last unary node. It's an optimization to prevent recursion in case that
			// we found two or more unary expressions after each other. For example the
			// expression "-!-1" contains only unary operators that will be parsed by
			// a single `parseExpression()` call.
			std::shared_ptr<AstUnary> lastUnaryNode = nullptr;
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

					std::string symbolName(_tokenizer.getTokenName(token));
					Token token_tmp = token;

					std::shared_ptr<AstNode> newNode;

					if (_tokenizer.next(&token) == TokenType::kTokenLParen)
					{
						// Will be parsed by `parseCall()` again.
						_tokenizer.set(&token_tmp);
						MATHPRESSO_PROPAGATE(parseCall(&newNode));
					}
					else
					{
						_tokenizer.set(&token);

						resolver::ContextPtr ctxfound;
						std::shared_ptr<AstSymbol> sym = resolver::resolveVariable(_shadowContext, symbolName, &ctxfound);

						if (sym)
						{
							if (!sym->isDeclared())
								MATHPRESSO_PARSER_ERROR(token, "Can't use variable '%s' that is being declared.", sym->getName());

							// Put symbol to shadow scope if it's global. This is done lazily and
							// only once per symbol when it's referenced.
							if (ctxfound->isGlobal())
							{
								sym = _ast->shadowSymbol(sym);
								MATHPRESSO_NULLCHECK(sym);

								sym->setVarSlotId(_ast->newSlotId());
							}

							newNode = std::make_shared<AstVar>();
							MATHPRESSO_NULLCHECK(newNode);
							std::static_pointer_cast<AstVar>(newNode)->setSymbol(sym);

							if (sym->hasSymbolFlag(AstSymbolFlags::kAstSymbolIsComplex))
								newNode->addNodeFlags(AstNodeFlags::kAstTakesComplex | AstNodeFlags::kAstReturnsComplex);

							newNode->setPosition(token.getPosAsUInt());
							sym->incUsedCount();
						}
						else
						{
							MATHPRESSO_PARSER_ERROR(token_tmp, "Unresolved symbol %s.", symbolName.c_str());
						}
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
					std::shared_ptr<AstImm> newNode = std::make_shared<AstImm>();
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
					std::shared_ptr<AstNode> newNode;
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
					std::string name(_tokenizer.getTokenName(token));
					auto op = resolver::resolveFunction(_shadowContext, name, 1);
					if (!op)
					{
						// for expressions like '----x'
						for (size_t i = 0; i < name.size(); i++)
						{
							if (op = resolver::resolveFunction(_shadowContext, name.substr(i, 1), 1))
							{
								std::shared_ptr<AstUnaryOp> opNode = std::make_shared<AstUnaryOp>();
								MATHPRESSO_NULLCHECK(opNode);
								opNode->setPosition(token.getPosAsUInt() + uint32_t(i));

								opNode->_mpOp = op;
								opNode->_opName = name.substr(i, 1);

								if (lastUnaryNode == nullptr)
									currentNode = opNode;
								else
									lastUnaryNode->setChild(opNode);

								isNested = true;
								lastUnaryNode = opNode;
							}
							else
							{
								MATHPRESSO_PARSER_ERROR(token, "Invalid unary operator %s.", name.c_str());
							}
						}

					}
					else
					{
						// Parse the unary operator.
						std::shared_ptr<AstUnaryOp> opNode = std::make_shared<AstUnaryOp>();
						MATHPRESSO_NULLCHECK(opNode);
						opNode->setPosition(token.getPosAsUInt());

						opNode->_mpOp = op;
						opNode->_opName = name;

						if (lastUnaryNode == nullptr)
							currentNode = opNode;
						else
							lastUnaryNode->setChild(opNode);

						isNested = true;
						lastUnaryNode = opNode;
					}
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
						currentNode = rootBinaryNode;
					}

					*pNodeOut = currentNode;
					return ErrorCode::kErrorOk;
				}

				// Parse Binary Operators
				case TokenType::kTokenOperator:
				{
					std::string name(_tokenizer.getTokenName(token));
					auto op = resolver::resolveFunction(_shadowContext, name, 2);

					if (name == "=")
					{
						op = resolver::resolveFunction(_shadowContext, "=", 1);
						// Check whether the assignment is valid.
						if (currentNode->getNodeType() != AstNodeType::kAstNodeVar)
							MATHPRESSO_PARSER_ERROR(token, "Can't assign to a non-variable.");

						std::shared_ptr<AstSymbol> sym = std::static_pointer_cast<AstVar>(currentNode)->getSymbol();
						if (sym->hasSymbolFlag(AstSymbolFlags::kAstSymbolIsReadOnly))
							MATHPRESSO_PARSER_ERROR(token, "Can't assign to a read-only variable '%s'.", sym->getName());

						if (isNested)
							MATHPRESSO_PARSER_ERROR(token, "Invalid assignment inside an expression.");

						sym->incWriteCount();
					}

					if (!op)
						MATHPRESSO_PARSER_ERROR(token, "Invalid Operator: " + name);

					std::shared_ptr<AstBinaryOp> newNode = std::make_shared<AstBinaryOp>();
					MATHPRESSO_NULLCHECK(newNode);

					newNode->_mpOp = op;
					newNode->_opName = name;

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
						rootBinaryNode = newNode;
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
							currentBinaryNode = std::static_pointer_cast<AstBinaryOp>(currentBinaryNode->getParent());
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
							rootBinaryNode = newNode;
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
							std::shared_ptr<AstNode> pNode = currentBinaryNode->unlinkRight();
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
	Error Parser::parseCall(std::shared_ptr<AstNode>* pNodeOut)
	{
		Token token;
		uint32_t uToken;

		uToken = _tokenizer.next(&token);
		MATHPRESSO_ASSERT(uToken == TokenType::kTokenSymbol);
		uint32_t position = token.getPosAsUInt();

		std::string fnName(_tokenizer.getTokenName(token));

		uToken = _tokenizer.next(&token);
		if (uToken != TokenType::kTokenLParen)
			MATHPRESSO_PARSER_ERROR(token, "Expected a '(' token after a function name.");

		std::shared_ptr<AstCall> callNode = std::make_shared<AstCall>();
		MATHPRESSO_NULLCHECK(callNode);

		callNode->setPosition(position);

		uToken = _tokenizer.peek(&token);
		if (uToken != TokenType::kTokenRParen)
		{
			for (;;)
			{
				// Parse the argument expression.
				std::shared_ptr<AstNode> expression;
				Error err = parseExpression(&expression, true);

				if (err != ErrorCode::kErrorOk)
				{
					return err;
				}

				callNode->appendNode(expression);

				// Parse ')' or ',' tokens.
				uToken = _tokenizer.peek(&token);
				if (uToken == TokenType::kTokenRParen)
					break;

				if (uToken != TokenType::kTokenComma)
				{
					MATHPRESSO_PARSER_ERROR(token, "Expected either ',' or ')' token.");
				}

				_tokenizer.consume();
			}
		}

		_tokenizer.consume();

		if (resolver::resolveFunction(_shadowContext, fnName, callNode->getLength()))
		{
			// the correct function to use is determined by the optimizer.
			callNode->_opName = fnName;
		}
		else
		{
			// TODO: object translation and error reporting possible signatures.
			MATHPRESSO_PARSER_ERROR(token, "Object '%s': wrong number of arguments. (Received %d).", fnName.c_str(), callNode->getLength());
		}

		*pNodeOut = callNode;
		return ErrorCode::kErrorOk;
	}

	Error Parser::reparseTernary(std::shared_ptr<AstNode> node)
	{
		for (size_t i = 0; i < node->getLength(); i++)
		{
			unsigned int ret;
			if (node->getNodeType() == AstNodeType::kAstNodeBinaryOp && node->_opName == "?")
			{

				std::shared_ptr<AstBinaryOp> lastColon = std::static_pointer_cast<AstBinaryOp>(node);
				// go to the last Colon after question-marks.
				while (lastColon->_opName == "?")
				{
					lastColon = std::static_pointer_cast<AstBinaryOp>(lastColon->getRight());
				}

				while (lastColon->getRight() && lastColon->getRight()->_opName == ":")
				{
					lastColon = std::static_pointer_cast<AstBinaryOp>(lastColon->getRight());
				}

				if (lastColon->_opName != ":")
				{
					return _errorReporter->onError(ErrorCode::kErrorInvalidSyntax, node->getPosition(),
												   "Invalid ternary operation. Expected a ':'.");
				}

				std::shared_ptr<AstNode> branchCondition = std::static_pointer_cast<AstBinaryOp>(node)->getLeft();
				std::shared_ptr<AstNode> branchLeft = lastColon->getLeft();
				std::shared_ptr<AstNode> branchRight = lastColon->getRight();

				// remove branchCondition from the AST
				std::static_pointer_cast<AstBinaryOp>(node)->setLeft(nullptr);
				branchCondition->_parent.reset();

				// remove the right path from the AST.
				lastColon->setRight(nullptr);
				branchRight->_parent.reset();


				// Distinguish between a complex and a non-complex case:
				// i.e.: cond1 ? cond2 ? a : b : c
				if (std::static_pointer_cast<AstBinaryOp>(node)->getRight() != lastColon)
				{
					// remove left branch from the AST.
					branchLeft = std::static_pointer_cast<AstBinaryOp>(node)->getRight();
					std::static_pointer_cast<AstBinaryOp>(node)->setRight(nullptr);
					branchLeft->_parent.reset();

					// correct the right path.
					std::shared_ptr<AstBinaryOp> preLastColon = std::static_pointer_cast<AstBinaryOp>(lastColon->getParent());
					preLastColon->replaceAt(1, lastColon->getLeft());
					lastColon->setLeft(nullptr);
				}
				// i.e.: cond1 ? a : b
				else
				{
					// remove left branch from the AST.
					lastColon->setLeft(nullptr);
					branchLeft->_parent.reset();
				}

				// create the new Ternary Node.
				std::shared_ptr<AstTernaryOp> ternaryNode = std::make_shared<AstTernaryOp>();
				ternaryNode->setCondition(branchCondition);
				ternaryNode->setLeft(branchLeft);
				ternaryNode->setRight(branchRight);
				ternaryNode->_mpOp = resolver::resolveFunction(_shadowContext, "_ternary_", 3);
				ternaryNode->_opName = "_ternary_";

				// add the new node to the AST.
				node->getParent()->replaceNode(node, ternaryNode);

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
					ret = ErrorCode::kErrorInvalidState;
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
