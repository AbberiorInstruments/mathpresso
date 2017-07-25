// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Export]
#define MATHPRESSO_EXPORTS

// [Dependencies]
#include <mathpresso/mpparser_p.h>

namespace mathpresso {

// ============================================================================
// [mathpresso::Parser - Syntax Error]
// ============================================================================

#define MATHPRESSO_PARSER_ERROR(_Token_, ...) \
  return _errorReporter->onError( \
    kErrorInvalidSyntax, static_cast<uint32_t>((size_t)(_Token_.position)), __VA_ARGS__)

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
    : AstScope(parser->_ast, parser->_currentScope, kAstScopeNested),
      _parser(parser) {
    _parser->_currentScope = this;
  }

  MATHPRESSO_INLINE ~AstNestedScope() {
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

Error Parser::parseProgram(AstProgram* block) {
  for (;;) {
    Token token;
    uint32_t uToken = _tokenizer.peek(&token);

    // Parse the end of the input.
    if (uToken == kTokenEnd)
      break;
    MATHPRESSO_PROPAGATE(parseStatement(block, kEnableVarDecls | kEnableNestedBlock));
  }

  if (block->getLength() == 0)
    return MATHPRESSO_TRACE_ERROR(kErrorNoExpression);

  return kErrorOk;
}

// Parse <statement>; or { [<statement>; ...] }
Error Parser::parseStatement(AstBlock* block, uint32_t flags) {
  Token token;
  uint32_t uToken = _tokenizer.peek(&token);

  // Parse the ';' token.
  if (uToken == kTokenSemicolon) {
    _tokenizer.consume();
    return kErrorOk;
  }
  
  if (uToken == kTokenColon) {
	  _tokenizer.consume();
	  return kErrorOk;
  }

  // Parse a nested block.
  if (uToken == kTokenLCurl) {
    AstBlock* nested;

    if (!(flags & kEnableNestedBlock))
      MATHPRESSO_PARSER_ERROR(token, "Cannot declare a new block-scope here.");

    MATHPRESSO_PROPAGATE(block->willAdd());
    MATHPRESSO_NULLCHECK(nested = _ast->newNode<AstBlock>());
    block->appendNode(nested);

    AstNestedScope tmpScope(this);
    return parseBlockOrStatement(nested);
  }

  // Parse a variable declaration.
  if (uToken == kTokenVar) {
    if (!(flags & kEnableVarDecls))
      MATHPRESSO_PARSER_ERROR(token, "Cannot declare a new variable here.");
    return parseVariableDecl(block);
  }

  // Parse an expression.
  AstNode* expression;

  MATHPRESSO_PROPAGATE(block->willAdd());
  MATHPRESSO_PROPAGATE(parseExpression(&expression, false));
  block->appendNode(expression);

  uToken = _tokenizer.peek(&token);
  if (uToken == kTokenSemicolon) {
    _tokenizer.consume();
    return kErrorOk;
  }

  if (uToken == kTokenColon) {
	  _tokenizer.consume();
	  return kErrorOk;
  }

  if (uToken == kTokenEnd)
    return kErrorOk;

  MATHPRESSO_PARSER_ERROR(token, "Expected a ';' after an expression.");
}

// Parse <block|statement>;.
Error Parser::parseBlockOrStatement(AstBlock* block) {
  Token token;
  uint32_t uToken = _tokenizer.next(&token);

  // Parse the <block>, consume '{' token.
  block->setPosition(token.getPosAsUInt());
  if (uToken == kTokenLCurl) {
    for (;;) {
      uToken = _tokenizer.peek(&token);

      // Parse the end of the block '}'.
      if (uToken == kTokenRCurl) {
        // Consume the '}' token.
        _tokenizer.consume();
        return kErrorOk;
      }
      else {
        MATHPRESSO_PROPAGATE(parseStatement(block, kEnableVarDecls | kEnableNestedBlock));
      }
    }
  }
  else {
    return parseStatement(block, kNoFlags);
  }
}

// Parse "var <name> = <expression>[, <name> = <expression>, ...];".
Error Parser::parseVariableDecl(AstBlock* block) {
  Token token;
  uint32_t uToken = _tokenizer.next(&token);
  StringRef str;

  bool isFirst = true;
  uint32_t position = token.getPosAsUInt();

  // Parse the 'var' keyword.
  if (uToken != kTokenVar)
    MATHPRESSO_PARSER_ERROR(token, "Expected 'var' keyword.");

  AstScope* scope = _currentScope;
  for (;;) {
    // Parse the variable name.
    if (_tokenizer.next(&token) != kTokenSymbol)
      MATHPRESSO_PARSER_ERROR(token, isFirst
        ? "Expected a variable name after 'var' keyword."
        : "Expected a variable name after colon ','.");

    MATHPRESSO_PROPAGATE(block->willAdd());

    if (!isFirst)
      position = token.getPosAsUInt();

    // Resolve the variable name.
    AstSymbol* vSym;
    AstScope* vScope;

    str.set(_tokenizer._start + token.position, token.length);
    if ((vSym = scope->resolveSymbol(str, token.hVal, &vScope)) != nullptr) {
      if (vSym->getSymbolType() != kAstSymbolVariable || scope == vScope)
        MATHPRESSO_PARSER_ERROR(token, "Attempt to redefine '%s'.", vSym->getName());

      if (vSym->hasNode()) {
        uint32_t line, column;
        _errorReporter->getLineAndColumn(vSym->getNode()->getPosition(), line, column);
        MATHPRESSO_PARSER_WARNING(token, "Variable '%s' shadows a variable declared at [%d:%d].", vSym->getName(), line, column);
      }
      else {
        MATHPRESSO_PARSER_WARNING(token, "Variable '%s' shadows a variable of the same name.", vSym->getName());
      }
    }

    vSym = _ast->newSymbol(str, token.hVal, kAstSymbolVariable, scope->getScopeType());
    MATHPRESSO_NULLCHECK(vSym);
    scope->putSymbol(vSym);

    AstVarDecl* decl = _ast->newNode<AstVarDecl>();
    MATHPRESSO_NULLCHECK_(decl, { _ast->deleteSymbol(vSym); });
	decl->mpOp_ = _ops->at("=$2").get();

    decl->setPosition(position);
    decl->setSymbol(vSym);

    // Assign a slot and fill to safe defaults.
    vSym->setVarOffset(0);
    vSym->setVarSlotId(_ast->newSlotId());
    vSym->setNode(decl);

    // Parse possible assignment '='.
    uToken = _tokenizer.next(&token);
    bool isAssigned = (uToken == kTokenAssign);

    if (isAssigned) {
      AstNode* expression;
      MATHPRESSO_PROPAGATE_(parseExpression(&expression, false), { _ast->deleteNode(decl); });
	
      decl->setChild(expression);
      vSym->incWriteCount();

      uToken = _tokenizer.next(&token);
    }

    // Make the symbol declared so it can be referenced after now.
    vSym->setDeclared();

    // Parse the ',' or ';' tokens.
    if (uToken == kTokenComma || uToken == kTokenSemicolon || uToken == kTokenEnd) {
      block->appendNode(decl);

      // Token ';' terminates the declaration.
      if (uToken != kTokenComma)
        break;
    }
    else {
      _ast->deleteSymbol(vSym);
      MATHPRESSO_PARSER_ERROR(token, "Unexpected token %d.", token.token);
    }

    isFirst = false;
  }

  return kErrorOk;
}

Error Parser::parseExpression(AstNode** pNode, bool isNested) {
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
  uint32_t op;

  // Current binary operator node. Initial nullptr value means that the parsing
  // just started and there is no binary operator yet. Once the first binary
  // operator has been parsed `currentBinaryNode` will be set accordingly.
  AstBinaryOp* currentBinaryNode = nullptr;

  // Currently parsed node.
  AstNode* currentNode = nullptr;
  
  for (;;) {
    // Last unary node. It's an optimization to prevent recursion in case that
    // we found two or more unary expressions after each other. For example the
    // expression "-!-1" contains only unary operators that will be parsed by
    // a single `parseExpression()` call.
    AstUnary* lastUnaryNode = nullptr;
	bool b_complex = false;

_Repeat1:
    switch (_tokenizer.next(&token)) {
      // Parse a variable, a constant, or a function-call. This can be repeated
      // one or several times based on the expression type. For unary nodes it's
      // repeated immediately, for binary nodes it's repeated after the binary
      // node is created.

      // Parse a symbol (variable or function name).
      case kTokenSymbol: {
        StringRef str(_tokenizer._start + token.position, token.length);

        AstScope* symScope;
        AstSymbol* sym = scope->resolveSymbol(str, token.hVal, &symScope);

        if (sym == nullptr)
          MATHPRESSO_PARSER_ERROR(token, "Unresolved symbol %.*s.", static_cast<int>(str.getLength()), str.getData());

        uint32_t symType = sym->getSymbolType();
        AstNode* newNode;

		if (symType == kAstSymbolVariable) {
			if (!sym->isDeclared())
				MATHPRESSO_PARSER_ERROR(token, "Can't use variable '%s' that is being declared.", sym->getName());

			// Put symbol to shadow scope if it's global. This is done lazily and
			// only once per symbol when it's referenced.
			if (symScope->isGlobal()) {
				sym = _ast->shadowSymbol(sym);
				MATHPRESSO_NULLCHECK(sym);

				sym->setVarSlotId(_ast->newSlotId());
				symScope = _ast->getRootScope();
				symScope->putSymbol(sym);
			}

			newNode = _ast->newNode<AstVar>();
			MATHPRESSO_NULLCHECK(newNode);
			static_cast<AstVar*>(newNode)->setSymbol(sym);

			if (sym->hasSymbolFlag(kAstSymbolIsComplex))
				newNode ->addNodeFlags(kAstTakesComplex | kAstReturnsComplex);

			newNode->setPosition(token.getPosAsUInt());
			sym->incUsedCount();
		}
        else {
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
	  case kTokenComplex:
		  b_complex = true;
      case kTokenNumber: {
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
      case kTokenComma:
      case kTokenSemicolon:
      case kTokenRParen: {
        MATHPRESSO_PARSER_ERROR(token, "Expected an expression.");
      }

      // Parse a nested expression.
      case kTokenLParen: {
        uint32_t position = token.getPosAsUInt();

        AstNode* newNode;
        MATHPRESSO_PROPAGATE(parseExpression(&newNode, true));

        if (_tokenizer.next(&token) != kTokenRParen)
          MATHPRESSO_PARSER_ERROR(token, "Expected a ')' token.");

        if (lastUnaryNode == nullptr)
          currentNode = newNode;
        else
          lastUnaryNode->setChild(newNode);

        break;
      }

      // Parse a right-to-left associative unary operator ('+', '-', "!").
      case kTokenAdd         : op = kOpNone  ; goto _Unary;
      case kTokenSub         : op = kOpNeg   ; goto _Unary;
      case kTokenNot         : op = kOpNot   ; goto _Unary;
_Unary: {

        // Parse the unary operator.
        AstUnaryOp* opNode = _ast->newNode<AstUnaryOp>(op);
        MATHPRESSO_NULLCHECK(opNode);
        opNode->setPosition(token.getPosAsUInt());

		std::string opNameDecorated(std::string(_tokenizer._start).substr(token.position, token.length) + "$1");

		if (_ops->find(opNameDecorated) != _ops->end())
		{
			opNode->mpOp_ = _ops->at(opNameDecorated).get();
		}

        if (lastUnaryNode == nullptr)
          currentNode = opNode;
        else
          lastUnaryNode->setChild(opNode);

        isNested = true;
        lastUnaryNode = opNode;

        goto _Repeat1;
      }

      case kTokenEnd: {
        MATHPRESSO_PARSER_ERROR(token, "Unexpected end of the program.");
      }

      default: {
        MATHPRESSO_PARSER_ERROR(token, "Unexpected token %d.", token.token);
      }
    }

// _Repeat2:
    switch (_tokenizer.next(&token)) {
      // Parse the expression terminators - ',', ';', ')' or EOI.
      case kTokenComma:
      case kTokenSemicolon:
      case kTokenRParen:
      case kTokenEnd: {
        _tokenizer.set(&token);

        if (currentBinaryNode != nullptr) {
          currentBinaryNode->setRight(currentNode);
          // Iterate to the top-most node.
          while (currentBinaryNode->hasParent())
            currentBinaryNode = static_cast<AstBinaryOp*>(currentBinaryNode->getParent());
          currentNode = currentBinaryNode;
        }

        *pNode = currentNode;
        return kErrorOk;
      }

      // Parse a binary operator.
      case kTokenAssign: {
        op = kOpAssign;
		
        // Check whether the assignment is valid.
        if (currentNode->getNodeType() != kAstNodeVar)
          MATHPRESSO_PARSER_ERROR(token, "Can't assign to a non-variable.");

        AstSymbol* sym = static_cast<AstVar*>(currentNode)->getSymbol();
        if (sym->hasSymbolFlag(kAstSymbolIsReadOnly))
          MATHPRESSO_PARSER_ERROR(token, "Can't assign to a read-only variable '%s'.", sym->getName());

        if (isNested)
          MATHPRESSO_PARSER_ERROR(token, "Invalid assignment inside an expression.");

        sym->incWriteCount();
        goto _Binary;
      }


      case kTokenEq          : op = kOpEq          ; goto _Binary;
      case kTokenNe          : op = kOpNe          ; goto _Binary;
      case kTokenGt          : op = kOpGt          ; goto _Binary;
      case kTokenGe          : op = kOpGe          ; goto _Binary;
      case kTokenLt          : op = kOpLt          ; goto _Binary;
      case kTokenLe          : op = kOpLe          ; goto _Binary;
      case kTokenAdd         : op = kOpAdd         ; goto _Binary;
      case kTokenSub         : op = kOpSub         ; goto _Binary;
      case kTokenMul         : op = kOpMul         ; goto _Binary;
      case kTokenDiv         : op = kOpDiv         ; goto _Binary;
      case kTokenMod         : op = kOpMod         ; goto _Binary;
	  case kTokenQMark       : op = kOpQMark       ; goto _Binary;
	  case kTokenColon       : op = kOpColon       ; goto _Binary;
_Binary: {
        AstBinaryOp* newNode = _ast->newNode<AstBinaryOp>(op);
        MATHPRESSO_NULLCHECK(newNode);

		std::string opNameDecorated(std::string(_tokenizer._start).substr(token.position, token.length) + "$2");
		
		if (_ops->find(opNameDecorated) != _ops->end())
		{
			newNode->mpOp_ = _ops->at(opNameDecorated).get();
		}

        newNode->setPosition(token.getPosAsUInt());

        if (currentBinaryNode == nullptr) {
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

        uint32_t currentBinaryPrec = OpInfo::get(currentBinaryNode->getOp()).precedence;
        uint32_t newBinaryPrec = OpInfo::get(op).precedence;

        if (currentBinaryPrec > newBinaryPrec) {
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
        else {
          currentBinaryNode->setRight(currentNode);

          // Advance to the top-most binaryNode that has less or equal precedence
          // than newBinaryPrec.
          while (currentBinaryNode->hasParent()) {
            // Terminate conditions:
            //   1. currentBinaryNode has higher precedence than newNode.
            //   2. currentBinaryNode has equal precedence and right-to-left associativity.
            if (OpInfo::get(currentBinaryNode->getOp()).rightAssociate(newBinaryPrec))
              break;
            currentBinaryNode = static_cast<AstBinaryOp*>(currentBinaryNode->getParent());
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
          if (!currentBinaryNode->hasParent() && !OpInfo::get(currentBinaryNode->getOp()).rightAssociate(newBinaryPrec)) {
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
          else {
            AstNode* pNode = currentBinaryNode->unlinkRight();
            currentBinaryNode->setRight(newNode);
            newNode->setLeft(pNode);
          }

          isNested = true;
          currentBinaryNode = newNode;

          break;
        }
      }

	 

      default: {
        MATHPRESSO_PARSER_ERROR(token, "Unexpected token %d.", token.token);
      }
    }
  }
}

// Parse "function([arg1 [, arg2, ...] ])".
Error Parser::parseCall(AstNode** pNodeOut) {
  Token token;
  uint32_t uToken;

  uToken = _tokenizer.next(&token);
  MATHPRESSO_ASSERT(uToken == kTokenSymbol);
  uint32_t position = token.getPosAsUInt();

  StringRef str(_tokenizer._start + token.position, token.length);
  AstSymbol* sym = _currentScope->resolveSymbol(str, token.hVal); // resolve the Symbol.

  if (sym == nullptr)
    MATHPRESSO_PARSER_ERROR(token, "Unresolved symbol %.*s.", static_cast<int>(str.getLength()), str.getData());

  if (sym->getSymbolType() != kAstSymbolIntrinsic &&
      sym->getSymbolType() != kAstSymbolFunction)
    MATHPRESSO_PARSER_ERROR(token, "Expected a function name.");

  uToken = _tokenizer.next(&token);
  if (uToken != kTokenLParen)
    MATHPRESSO_PARSER_ERROR(token, "Expected a '(' token after a function name.");

  AstCall* callNode = _ast->newNode<AstCall>();
  MATHPRESSO_NULLCHECK(callNode);

  callNode->setSymbol(sym); // set symbol as part of the node.
  callNode->setPosition(position);

  uToken = _tokenizer.peek(&token);
  if (uToken != kTokenRParen) { // append parameters as children (sym not necessary)
    for (;;) {
      // Parse the argument expression.
      AstNode* expression;
      Error err;

      if ((err = callNode->willAdd()) != kErrorOk || (err = parseExpression(&expression, true)) != kErrorOk) {
        _ast->deleteNode(callNode);
        return err;
      }

      callNode->appendNode(expression);

      // Parse ')' or ',' tokens.
      uToken = _tokenizer.peek(&token);
      if (uToken == kTokenRParen)
        break;

      if (uToken != kTokenComma) {
        _ast->deleteNode(callNode);
        MATHPRESSO_PARSER_ERROR(token, "Expected either ',' or ')' token.");
      }

      _tokenizer.consume();
    }
  }

  _tokenizer.consume();

  // Validate the number of function arguments.
  size_t n = callNode->getLength(); 
  // lookup symbol here? _crrentScope wrong value? 
  uint32_t reqArgs = sym->getFuncArgs();

  if (n != reqArgs) {
    _ast->deleteNode(callNode);
    MATHPRESSO_PARSER_ERROR(token, "Function '%s' requires %u argument(s) (%u provided).", sym->getName(), reqArgs, n);
  }

  std::string opNameDecorated(str.getData(), str.getLength());
  opNameDecorated += "$" + std::to_string(n);

  if (_ops->find(opNameDecorated) != _ops->end())
  {
	  callNode->_symbol->setOp(_ops->at(opNameDecorated));
  }
  else
  {
	  MATHPRESSO_PARSER_ERROR(token, "Function '%s' requires a MpOperation-Object, which is not to be found.", sym->getName());
  }

  // Transform an intrinsic function into unary or binary operator.
  if (sym->getSymbolType() == kAstSymbolIntrinsic) {
    const OpInfo& op = OpInfo::get(sym->getOpType());
    MATHPRESSO_ASSERT(n == op.getOpCount());

    AstNode* opNode;
    if (reqArgs == 1) {
      AstUnaryOp* unary = _ast->newNode<AstUnaryOp>(op.type);
      MATHPRESSO_NULLCHECK(unary);

	  unary->mpOp_ = callNode->getSymbol()->getOp().get();
	  
      unary->setChild(callNode->removeAt(0));
      opNode = unary;
    }
    else {
      AstBinaryOp* binary = _ast->newNode<AstBinaryOp>(op.type);
      MATHPRESSO_NULLCHECK(binary);

	  binary->mpOp_ = callNode->getSymbol()->getOp().get();

      binary->setRight(callNode->removeAt(1));
      binary->setLeft(callNode->removeAt(0));
      opNode = binary;
    }

    opNode->setPosition(callNode->getPosition());
    _ast->deleteNode(callNode);

    *pNodeOut = opNode;
    return kErrorOk;
  }
  else {
    *pNodeOut = callNode;
    return kErrorOk;
  }
}

} // mathpresso namespace
