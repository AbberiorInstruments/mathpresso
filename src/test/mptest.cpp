// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

#include "../mathpresso/mathpresso.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct TestExpression {
  const char* expression;
  double expected;
};

#define TABLE_SIZE(table) \
  (sizeof(table) / sizeof(table[0]))

#define TEST_EXPRESSION(expression) \
  { #expression, (double)(expression) }

#if defined(min)
 #undef min
#endif // min

#if defined(max)
# undef max
#endif // max

namespace environment {
  const double E  = 2.7182818284590452354;
  const double PI = 3.14159265358979323846;

  inline double avg(double x, double y) { return (x + y) * 0.5; }
  inline double min(double x, double y) { return x < y ? x : y; }
  inline double max(double x, double y) { return x > y ? x : y; }

  inline double recip(double x) { return 1.0 / x; }
}

int main(int argc, char* argv[]) {
  using namespace environment;

  bool failed = false;

  // Made volatile so the compiler doesn't optimize evaluation of the
  // expression into 80-bit FPU registers. Otherwise we will have some
  // failures.
  volatile double x = 5.1;
  volatile double y = 6.7;
  volatile double z = 9.9;

  mathpresso::Context ctx;
  mathpresso::Expression e0;
  mathpresso::Expression e1;

  ctx.addEnvironment(mathpresso::kMPEnvironmentAll);
  ctx.addVariable("x", 0 * sizeof(double));
  ctx.addVariable("y", 1 * sizeof(double));
  ctx.addVariable("z", 2 * sizeof(double));

  double variables[] = { x, y, z };

  TestExpression tests[] = {
    TEST_EXPRESSION( x+y ),
    TEST_EXPRESSION( x-y ),
    TEST_EXPRESSION( x*y ),
    TEST_EXPRESSION( x/y ),
    TEST_EXPRESSION( -(x+y) ),
    TEST_EXPRESSION( -(x-y) ),
    TEST_EXPRESSION( x*z + y*z),
    TEST_EXPRESSION( x*z - y*z),
    TEST_EXPRESSION( x*z*y*z),
    TEST_EXPRESSION( x*z/y*z),
    TEST_EXPRESSION( x == y ),
    TEST_EXPRESSION( x != y ),
    TEST_EXPRESSION( x <  y ),
    TEST_EXPRESSION( x <= y ),
    TEST_EXPRESSION( x >  y ),
    TEST_EXPRESSION( x >= y ),
    TEST_EXPRESSION( x + y == y - z ),
    TEST_EXPRESSION( x > y == y < z ),
    TEST_EXPRESSION( -x ),
    TEST_EXPRESSION( -1 + x ),
    TEST_EXPRESSION( -(-(-1)) ),
    TEST_EXPRESSION( -(-(-x)) ),
    TEST_EXPRESSION( (x+y)*x ),
    TEST_EXPRESSION( (x+y)*y ),
    TEST_EXPRESSION( (x+y)*(1.19+z) ),
    TEST_EXPRESSION( ((x+(x+2.13))*y) ),
    TEST_EXPRESSION( (x+y+z*2+(x*z+z*1.5)) ),
    TEST_EXPRESSION( (((((((x-0.28)+y)+x)+x)*x)/1.12)*y) ),
    TEST_EXPRESSION( ((((x*((((y-1.50)+1.82)-x)/PI))/x)*x)+z) ),
    TEST_EXPRESSION( (((((((((x+1.35)+PI)/PI)-y)+z)-z)+y)/x)+0.81) ),
    TEST_EXPRESSION( round(x) ),
    TEST_EXPRESSION( round(-x) ),
    TEST_EXPRESSION( floor(x) ),
    TEST_EXPRESSION( floor(-x) ),
    TEST_EXPRESSION( ceil(x) ),
    TEST_EXPRESSION( ceil(-x) ),
    TEST_EXPRESSION( abs(-x) ),
    TEST_EXPRESSION( sqrt(x) ),
    TEST_EXPRESSION( recip(x) ),
    TEST_EXPRESSION( exp(x) ),
    TEST_EXPRESSION( log(x) ),
    TEST_EXPRESSION( log10(x) ),
    TEST_EXPRESSION( sin(x) ),
    TEST_EXPRESSION( cos(x) ),
    TEST_EXPRESSION( tan(x) ),
    TEST_EXPRESSION( sin(x) * cos(y) * tan(z) ),
    TEST_EXPRESSION( avg(x, y) ),
    TEST_EXPRESSION( min(x, y) ),
    TEST_EXPRESSION( max(x, y) ),
    TEST_EXPRESSION( pow(x, y) )
  };

  for (int i = 0; i < TABLE_SIZE(tests); i++) {
    const char* exp = tests[i].expression;
    bool expOk = true;

    int err;
    if ((err = e0.create(ctx, exp, mathpresso::kMPOptionDisableJIT))) {
      printf("[Failure]: \"%s\" (eval)\n", exp);
      printf("           ERROR %d (jit disabled).\n", err);
      continue;
    }

    if ((err = e1.create(ctx, exp, mathpresso::kMPOptionNone))) {
      printf("[Failure]: \"%s\" (jit)\n", exp);
      printf("           ERROR %d (jit enabled).\n", err);
      continue;
    }

    double expected = tests[i].expected;
    double res0 = e0.evaluate(variables);
    double res1 = e1.evaluate(variables);

    if (res0 != expected) {
      printf("[Failure]: \"%s\" (eval)\n", exp);
      printf("           result(%f) != expected(%f)\n", res0, expected);
      expOk = false;
    }

    if (res1 != expected) {
      printf("[Failure]: \"%s\" (jit)\n", exp);
      printf("           result(%f) != expected(%f)\n", res1, expected);
      expOk = false;
    }

    if (expOk) {
      printf("[Success]: \"%s\" -> %f\n", exp, expected);
    }
    else {
      failed = true;
    }
  }

  return failed ? 1 : 0;
}
