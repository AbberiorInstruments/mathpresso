MathPresso
==========

Mathematical Expression Parser And JIT Compiler.

  * [Official Repository (kobalicek/mathpresso)](https://github.com/kobalicek/mathpresso)
  * [Official Blog (asmbits)] (https://asmbits.blogspot.com/ncr)
  * [Official Chat (gitter)](https://gitter.im/kobalicek/mpsl)
  * [Permissive ZLIB license](./LICENSE.md)

Introduction
------------
  
MathPresso is a C++ library designed to parse mathematical expressions and compile them into machine code. It's much faster than traditional AST or byte-code based evaluators, because there is basically no overhead in the expression's execution. The JIT compiler is based on AsmJit and works on X86 and X64 architectures.

This is an updated version of MathPresso that uses a stripped-off MPSL engine designed to work with scalar double precision floating points. It has many bugs fixed compared to the last version on google-code and contains improvements that can make execution of certain built-in functions (intrinsics) faster if the host CPU supports SSE4.1 (rounding, fraction, modulo, etc...).

This is also a transitional version that is available to users that want to use MathPresso and cannot wait for the new MPSL engine, which is a work in progress.

### Additions by Abberior:

We added support for complex numbers. The compiler can create functions that take and return complex values and calculate them as needed. Every Function in Features, that suports complex numbers is available as a real and a complex version, where the Optimizer chooses the correct one.

We allso added the posibility for named sub-contexts similar to namespaces in c++.

Features
--------

  * Unary operators:
    * Negate `-(x)`
    * Not `!(x)`
  * Arithmetic operators:
    * Assignment `x = y`
    * Addition `x + y`
    * Subtraction `x - y`
    * Multiplication `x * y`
    * Division `x / y`
    * Modulo `x % y`
  * Comparison operators:
    * Equal `x == y`
    * Not equal `x != y`
    * Greater `x > y`
    * Greater or equal `x >= y`
    * Lesser `x < y`
    * Lesser or equal `x <= y`
  * Functions:
    * Check for NaN `isnan(x)`
    * Check for infinity `isinf(x)`
    * Check for finite number `isfinite(x)`
    * Get a sign bit `signbit(x)`
    * Copy sign `copysign(x, y)`
    * Round to nearest `round(x)`
    * Round to even `roundeven(x)`
    * Truncate `trunc(x)`
    * Floor `floor(x)`
    * Ceil `ceil(x)`
    * Average value `avg(x, y)`
    * Minimum value `min(x, y)`
    * Maximum value `max(x, y)`
    * Absolute value `abs(x)`
    * Exponential `exp(x)`
    * Logarithm `log(x)`
    * Logarithm of base 2 `log2(x)`
    * Logarithm of base 10 `log10(x)`
    * Square root `sqrt(x)`
    * Square root (complex return) `sqrtc(x)`
    * Fraction `frac(x)`
    * Reciprocal `recip(x)`
    * Power `pow(x, y)`
    * Sine `sin(x)`
    * Cosine `cos(x)`
    * Tangent `tan(x)`
    * Hyperbolic sine `sinh(x)`
    * Hyperbolic cosine `cosh(x)`
    * Hyperbolic tangent `tanh(x)`
    * Arcsine `asin(x)`
    * Arccosine `acos(x)`
    * Arctangent `atan(x)`
    * Arctangent `atan2(x, y)`
    * Conjugation `conjug(x)`
    * Get real part of a complex value `real(x)`
    * Get imaginary part of a complex value `imag(x)`
  * Other:
    * Ternary Operator `<bool> ? <expression> : <expression>`
    * Declaration of variables `var <name> = <val>`
  * Constants defined by `addBuiltIns()`:
    * Infinity `INF`
    * Not a Number `NaN`
    * Euler's constant `E = 2.7182818284590452354`
    * PI `PI = 3.14159265358979323846`
    * Imaginary number `i = sqrt(-1)`

Usage
-----

MathPresso's expression is always created around a `mathpresso::Context`, which defines an environment the expression can access and use. For example if you plan to extend MathPresso with your own function or constant the `Context` is the way to go. The `Context` also defines inputs and outputs of the expression as shown in the example below:

```c++
#include <mathpresso/mathpresso.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
  std::shared_ptr<mathpresso::Context> ctx = std::make_shared<mathpresso::Context>();
  std::shared_ptr<mathpresso::Expression> exp = std::make_shared<mathpresso::Expression>();

  // Initialize the context by adding MathPresso built-ins. Without this line
  // the context has no Operations it can do.
  ctx->addBuiltIns();

  // Let the context know the name of the variables we will refer to and
  // their positions in the data pointer. We will use an array of 3 doubles,
  // so index them by using `sizeof(double)`, like a normal C array.
  //
  // The `addVariable()` also contains a third parameter that describes
  // variable flags, use `kVariableRO` to make a certain variable read-only.
  ctx->addVariable("x", 0 * sizeof(double));
  ctx->addVariable("y", 1 * sizeof(double));
  ctx->addVariable("z", 2 * sizeof(double));

  // Compile the expression.
  //
  // The create parameters are:
  //   1. `std::shared_ptr<mathpresso::context>` - The expression's context / environment.
  //   2. `std::string body` - The expression body.
  //   3. `unsigned int` - Options, just pass `mathpresso::kNoOptions`.
  mathpresso::Error err = exp->compile(ctx, "(x*y) % z", mathpresso::kNoOptions);

  // Handle possible syntax or compilation error.
  if (err != mathpresso::kErrorOk) {
    printf("Expression Error: %u\n", err);
    return 1;
  }

  // To evaluate the expression you need to create the `data` to be passed
  // to the expression and initialize it. Every expression returns `double`,
  // to return more members simply use the passed `data`.
  double data[] = {
    1.2, // 'x' - available at data[0]
    3.8, // 'y' - available at data[1]
    1.3  // 'z' - available at data[2]
  };
  printf("Output: %f\n", exp->evaluate(data));

  return 0;
}
```

The example above should be self-explanatory. The next example does the same but by using a `struct` instead of an array to address the expression's data:

```c++
#include <mathpresso/mathpresso.h>
#include <stdio.h>

struct Data {
  inline Data(double x, double y, double z)
    : x(x), y(y), z(z) {}

  double x, y, z;
};

int main(int argc, char* argv[]) {
  std::shared_ptr<mathpresso::Context> ctx = std::make_shared<mathpresso::Context>();
  std::shared_ptr<mathpresso::Expression> exp = std::make_shared<mathpresso::Expression>();

  ctx->addBuiltIns();
  ctx->addVariable("x", MATHPRESSO_OFFSET(Data, x));
  ctx->addVariable("y", MATHPRESSO_OFFSET(Data, y));
  ctx->addVariable("z", MATHPRESSO_OFFSET(Data, z));

  mathpresso::Error err = exp->compile(ctx, "(x*y) % z", mathpresso::kNoOptions);
  if (err != mathpresso::kErrorOk) {
    printf("Expression Error: %u\n", err);
    return 1;
  }

  Data data(1.2, 3.8. 1.3);
  printf("Output: %f\n", exp->evaluate(&data));

  return 0;
}
```

Creation of custom operations
-----------------------------

We created a way to add your own functions or operators to mathpresso. You should find the necessary interface to implement in `MpOperatin.h`.

Every operation in mathpresso is implementing this interface.

```c++
#include <mathpresso/mathpresso.h>
#include <mathpresso/MpOperation.h>
#include <stdio.h>


// An example of a function computing the absolute value of a real number.
class MpOperationAbs : public MpOperationEval<double, double>
{
public:
  MpOperationAbs() noexcept : MpOperationEval<double, double>(1)
  {
  }

  // The function that is computed, if the parameter is an 
  // imediate value, that is known at compiletime
  virtual double evaluate(const double * args) const override;

  // The function used to generate the machine-code for this function.
  // This is called by the JitCompiler.
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

int main(int argc, char* argv[]) {
  std::shared_ptr<mathpresso::Context> ctx = std::make_shared<mathpresso::Context>();

  // This way you can add a function from c++ to the context.
  ctx->addObject("sin",  _OBJ(static_cast<double(*)(double, double)>(std::sin));

  //Adding the Function we declared earlier:
  ctx->addObject("abs",  std::make_shared<MpOperationAbs>());

  ctx->addVariable("x", 0 * sizeof(double));

  mathpresso::Error err = exp->compile(ctx, "sin(abs(x))", mathpresso::kNoOptions);
  if (err != mathpresso::kErrorOk) {
    printf("Expression Error: %u\n", err);
    return 1;
  }

  double[] data = { 3.5 };
  printf("Output: %f\n", exp->evaluate(&data));

  return 0;
}
```


Hirarchical named context
-------------------------

Hirarchical contexts allow the separation of functions and variables into structures like namespaces from c++.

```c++
#include <mathpresso/mathpresso.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
  std::shared_ptr<mathpresso::Context> ctx = std::make_shared<mathpresso::Context>();
  std::shared_ptr<mathpresso::Context> subctx = std::make_shared<mathpresso::Context>();

  // We add subctx as `sub` to our cotext.
  ctx->addChild("sub", subctx);

  // Add functions to the contexts
  ctx->addObject("sin",  _OBJ(static_cast<double(*)(double, double)>(std::sin));
  subctx->addObject("cos",  _OBJ(static_cast<double(*)(double, double)>(std::cos));

  ctx->addVariable("x", 0 * sizeof(double));

  mathpresso::Error err = exp->compile(ctx, "sin(sub.cos(x))", mathpresso::kNoOptions);
  if (err != mathpresso::kErrorOk) {
    printf("Expression Error: %u\n", err);
    return 1;
  }

  double[] data = { 3.5 };
  printf("Output: %f\n", exp->evaluate(&data));

  return 0;
}
```


Error Handling
--------------

MathPresso allows to attach an `OutputLog` instance to retrieve a human readable error message in case of error. It can output the following:

  - Errors, only one as MathPresso stops after the first error
  - Warnings
  - Abstract syntax tree (AST)
  - Assembly (ASM)

Here is the minimum working example that uses `OutputLog` to display errors. The interface is very simple, but extensible.

```c++
// This is a minimum working example that uses most of MathPresso features. It
// shows how to compile and evaluate expressions and how to handle errors. It
// also shows how to print the generated AST and machine code.
#include <mathpresso/mathpresso.h>
#include <stdio.h>

// The data passed to the expression.
struct Data {
  double x, y, z;
};

// By inheriting `OutputLog` one can create a way how to handle possible errors
// and report them to humans. The most interesting and used message type is
// `kMessageError`, because it signalizes an invalid expression. Other message
// types are used mostly for debugging.
struct MyOutputLog : public mathpresso::OutputLog {
  MyOutputLog() {}
  virtual ~MyOutputLog() {}
  virtual void log(unsigned int type, unsigned int line, unsigned int column, const char* message, size_t len) {
    switch (type) {
      case kMessageError:
        printf("[ERROR]: %s (line %u, column %u)\n", message, line, column);
        break;

      case kMessageWarning:
        printf("[WARNING]: %s (line %u, column %u)\n", message, line, column);
        break;

      case kMessageAstInitial:
        printf("[AST-INITIAL]\n%s", message);
        break;

      case kMessageAstFinal:
        printf("[AST-FINAL]\n%s", message);
        break;

      case kMessageAsm:
        printf("[ASSEMBLY]\n%s", message);
        break;

      default:
        printf("[UNKNOWN]\n%s", message);
        break;
    }
  }
};

int main(int argc, char* argv[]) {
  MyOutputLog outputLog;

  // Create the context, add builtins and define the `Data` layout.
  std::shared_ptr<mathpresso::Context> ctx = std::make_shared<mathpresso::Context>();
  ctx->addBuiltIns();
  ctx->addVariable("x" , MATHPRESSO_OFFSET(Data, x));
  ctx->addVariable("y" , MATHPRESSO_OFFSET(Data, y));
  ctx->addVariable("z" , MATHPRESSO_OFFSET(Data, z));

  // The following options will cause that MathPresso will send everything
  // it does to `OutputLog`.
  unsigned int options = 
    mathpresso::kOptionVerbose  | // Enable warnings, not just errors.
    mathpresso::kOptionDebugAst | // Enable AST dumps.
    mathpresso::kOptionDebugAsm ; // Enable ASM dumps.

  std::shared_ptr<mathpresso::Expression> exp = std::make_shared<mathpresso::Expression>();
  mathpresso::Error err = exp->compile(ctx,
    "-(-(abs(x * y - floor(x)))) * z * (12.9 - 3)", options, &outputLog);

  // Handle possible syntax or compilation error. The OutputLog has already
  // received and printed the reason in a human readable form.
  if (err) {
    printf("ERROR %u\n", err);
    return 1;
  }

  // Evaluate the expression, if compiled.
  Data data = { 12.2, 9.2, -1.9 };
  double result = exp->evaluate(&data);

  printf("RESULT: %f\n", result);
  return 0;
}
```

When executed the output of the application would be something like:

```
[AST-INITIAL]
* [Binary, <real> -> <real>]
  * [Binary, <real> -> <real>]
    - [Unary, <real> -> <real>]
      - [Unary, <real> -> <real>]
       abs(), <real> -> <real>
         - [Binary, <real> -> <real>]
           * [Binary, <real> -> <real>]
             x <real>
             y <real>
           floor(), <real> -> <real>
             x <real>
    z <real>
  - [Binary, <real> -> <real>]
    12.900000, <real>
    3.000000, <real>

[AST-FINAL]
* [Binary, <real> -> <real>]
  * [Binary, <real> -> <real>]
    abs(), <real> -> <real>
      - [Binary, <real> -> <real>]
        * [Binary, <real> -> <real>]
          x <real>
          y <real>
        floor(), <real> -> <real>
          x <real>
    z <real>
  9.900000, <real>

[ASSEMBLY]
L0:
lea rax, [L2]                           ; 488D05........          | lea pConst, [L2]
movsd xmm0, [rdx]                       ; F20F1002                | movsd v3, [pVariables]
mulsd xmm0, [rdx+8]                     ; F20F594208              | mulsd v3, [pVariables+8]
movsd xmm1, [rdx]                       ; F20F100A                | movsd v4, [pVariables]
roundsd xmm2, xmm1, 9                   ; 660F3A0BD109            | roundsd v5, v4, 9
subsd xmm0, xmm2                        ; F20F5CC2                | subsd v3, v5
xorpd xmm1, xmm1                        ; 660F57C9                | xorpd v6, v6
subsd xmm1, xmm0                        ; F20F5CC8                | subsd v6, v3
maxsd xmm1, xmm0                        ; F20F5FC8                | maxsd v6, v3
mulsd xmm1, [rdx+16]                    ; F20F594A10              | mulsd v6, [pVariables+16]
mulsd xmm1, [rax]                       ; F20F5908                | mulsd v6, [pConst]
movsd [rcx], xmm1                       ; F20F1109                | movsd [pResult], v6
L1:
ret                                     ; C3
.align 8
L2:
.data CDCCCCCCCCCC2340

RESULT: -1885.514400
```

Dependencies
------------

  * AsmJit - 1.0 or later.

Authors & Maintainers
---------------------

  * Petr Kobalicek <kobalicek.petr@gmail.com>
  * for Abberior Instruments.
    * Simon Leidenbach
