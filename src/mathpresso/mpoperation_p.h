// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MP_OPERATION_P_H
#define _MP_OPERATION_P_H

#include <complex>
#include "mpast_p.h"
#include "mpcompiler_p.h"
#include "mpoptimizer_p.h"

namespace mathpresso {

	typedef JitVar(*mpAsmFunc)(JitCompiler*, JitVar*);

	enum OperationFlags 
	{
		OpFlagNone = 0,
		OpFlagisOperator = 0x00000001,

		// Types of Operation that are allowed.
		OpFlagHasCtoC = 0x0000010,
		OpFlagHasCtoD = 0x0000020,
		OpFlagHasDtoC = 0x0000040,
		OpFlagHasDtoD = 0x0000080,
		OpFlagTypeMask = 0x00000F0,

		OpFlagHasCtoCasm = 0x0000100,
		OpFlagHasCtoDasm = 0x0000200,
		OpFlagHasDtoCasm = 0x0000400,
		OpFlagHasDtoDasm = 0x0000800,
		OpFlagTypeAsmMask = 0x0000F00,

		// binary Flags:
		OpFlagNopIfLZero = 0x10000000,
		OpFlagNopIfRZero = 0x20000000,
		OpFlagNopIfLOne = 0x40000000,
		OpFlagNopIfROne = 0x80000000,

		OpFlagNopIfZero = OpFlagNopIfLZero | OpFlagNopIfRZero,
		OpFlagNopIfOne = OpFlagNopIfLOne  | OpFlagNopIfROne
	};

	class OperationGeneric {
	public:
		// Con-/Destructor
		OperationGeneric(std::string name) : 
			numArgs_(0),
			flags_(OperationFlags::OpFlagNone),
			name_(name),
			fnCToC(nullptr),
			fnCToD(nullptr),
			fnDToC(nullptr),
			fnDToD(nullptr)
		{}

		virtual ~OperationGeneric() {}

		// Compile expression return true, if no error occurs, defaults to function calls.
		virtual JitVar compileToASM(JitCompiler* jc, AstNode * node);

		// Optimize in a more general way than with evaluate, i.e.: x * 1 -> x;
		// Right now: only flags are set and immediates are resolved.
		virtual Error optimize(AstOptimizer *opt, AstNode *node);

		// Evaluate, i.e. for Optimization of immediate values.
		virtual std::complex<double> evaluateCtoC(std::complex<double> *args);
		virtual std::complex<double> evaluateDtoC(double *args);
		virtual double evaluateCtoD(std::complex<double> *args);
		virtual double evaluateDtoD(double *args);

		uint32_t flags() { return flags_; }

	private:
		// Members:
		size_t numArgs_;
		uint32_t flags_;
		std::string name_;

		// Function-pointer:
		void * fnCToC;
		void * fnDToC;
		void * fnCToD;
		void * fnDToD;

		// Assembler-generators, eventually only in child classes.
		mpAsmFunc asmCToC;
		mpAsmFunc asmDToC;
		mpAsmFunc asmCToD;
		mpAsmFunc asmDToD;
	};
}


#endif //_MP_OPERATION_P_H