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
		OpFlagIsOperator = 0x00000001,
		OpFlagHasState = 0x00000002,

		// Types of Operation that are allowed.
		OpFlagCReturnsD = 0x0000010,
		OpFlagDReturnsC = 0x0000020,
		
		// binary Flags:
		OpFlagNopIfLZero = 0x10000000,
		OpFlagNopIfRZero = 0x20000000,
		OpFlagNopIfLOne  = 0x40000000,
		OpFlagNopIfROne  = 0x80000000,

		OpFlagNopIfZero = OpFlagNopIfLZero | OpFlagNopIfRZero,
		OpFlagNopIfOne  = OpFlagNopIfLOne  | OpFlagNopIfROne
	};

	class MpOperation 
	{
	public:
		// Con-/Destructor
		MpOperation(uint32_t nargs, uint32_t flags) :
			nargs_(nargs),
			flags_(flags)
		{}

		virtual ~MpOperation() 
		{
		}

		// Add ASM code to compiler stack 
		virtual JitVar compile(JitCompiler *jc, AstNode * node)  = 0;
		// Optimize AST 
		virtual Error optimize(AstOptimizer *opt, AstNode *node) = 0;

		uint32_t flags() 
		{ 
			return flags_; 
		}
		uint32_t nargs() 
		{ 
			return nargs_; 
		}

		void addFlags(uint32_t flags)
		{
			flags_ |= flags;
		}

	protected:
		uint32_t nargs_;
		uint32_t flags_;
	};

	class MpOperationFunc : public MpOperation
	{
	public:
		// Con-/Destructor
		MpOperationFunc(uint32_t nargs, uint32_t flags, void * fnD, void * fnC) : MpOperation(nargs, flags),
			fnD_(fnD),
			fnC_(fnC)
		{
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
		virtual Error optimize(AstOptimizer *opt, AstNode *node) override;
		virtual void setFn(void * fn, bool isComplex = false);
	protected:
		virtual double evaluateDRetD(double *args);
		// Function-pointer:
		void * fnC_;
		void * fnD_;
	};

	class MpOperationOp : public MpOperationFunc
	{
	public:
		MpOperationOp(uint32_t nargs, uint32_t flags, void * fnD, void * fnC, mpAsmFunc asmC, mpAsmFunc asmD) :
			MpOperationFunc(nargs, flags | OperationFlags::OpFlagIsOperator, fnD, fnC),
			asmC_(asmC),
			asmD_(asmD)
		{
		}
		
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
		virtual Error optimize(AstOptimizer *opt, AstNode *node) override;

	protected:
		mpAsmFunc asmC_;
		mpAsmFunc asmD_;
	};
}


#endif //_MP_OPERATION_P_H