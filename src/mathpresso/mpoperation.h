#pragma once
// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MP_OPERATION_P_H
#define _MP_OPERATION_P_H

#include <mathpresso/mathpresso.h>
#include <complex>
#include <limits>

namespace mathpresso {

	// Forward Declarations:
	struct JitCompiler;
	struct AstNode;
	struct AstOptimizer;
	struct JitVar;
	struct Context;

	typedef JitVar(*mpAsmFunc)(JitCompiler*, JitVar*);

	enum MpOperationFlags
	{
		OpFlagNone = 0,
		OpFlagIsOperator = 0x00000001,
		OpFlagHasState = 0x00000002,
		OpFlagHasAsm = 0x00000004,
		OpIsRighttoLeft = 0x00000008,

		OpIsCommutativ = 0x00000010,

		// Set, if no (complex|real) function is available.
		OpHasNoComplex = 0x00000020,
		OpHasNoReal = 0x00000040,

		// Types of Operation that are allowed.
		OpFlagCReturnsD = 0x0000100,
		OpFlagDReturnsC = 0x0000200,

		// some information, that is necessary for parsing:
		OpIsAssgignment = 0x0001000,

		// Set for optimization of binary operations
		OpFlagNopIfLZero = 0x10000000,
		OpFlagNopIfRZero = 0x20000000,
		OpFlagNopIfLOne = 0x40000000,
		OpFlagNopIfROne = 0x80000000,
		OpFlagNopIfZero = OpFlagNopIfLZero | OpFlagNopIfRZero,
		OpFlagNopIfOne  = OpFlagNopIfLOne | OpFlagNopIfROne
	};

	class MATHPRESSO_API MpOperation
	{
	public:
		struct Signature
		{
			enum class type
			{
				real = 0,
				complex = 1
			};

			Signature(size_t nargs) :
				return_type_(type::real),
				parameters_(nargs, { type::real, "" })
			{
			}

			type return_type_;

			struct param
			{
				type type_;
				std::string name_;
			};

			std::vector<param> parameters_;
		};


		// Con-/Destructor
		MpOperation(const Signature &s, uint32_t priority = 0) :
			priority_(priority)
		{
			nargs_ = uint32_t(s.parameters_.size());
			// Todo: set flags_
		}

		virtual ~MpOperation()
		{
		}

		// Add ASM code to compiler stack 
		virtual JitVar compile(JitCompiler *jc, AstNode * node) = 0;
		// Optimize AST 
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) = 0;

		uint32_t nargs() const 
		{ 
			return nargs_; 
		}
		bool isRightToLeft() const 
		{ 
			return (flags_ & OpIsRighttoLeft) != 0; 
		}
		uint32_t precedence() const 
		{ 
			return priority_; 
		}
		uint32_t flags() const 
		{ 
			return flags_; 
		}
	protected:
		uint32_t nargs_;
		uint32_t flags_;
		uint32_t priority_;
	};

	class MATHPRESSO_API MpOperationFunc : public MpOperation
	{
	public:
		// Con-/Destructor
		MpOperationFunc(uint32_t nargs, uint32_t flags, void * fnD, void * fnC) : MpOperation(nargs, flags),
			fnD_(fnD),
			fnC_(fnC)
		{
		}

		virtual ~MpOperationFunc()
		{}

		bool hasFlag(uint32_t flag)
		{
			return flag & flags_;
		}

		void addFlags(uint32_t flags)
		{
			flags_ |= flags;
		}

		void removeFlags(uint32_t flags)
		{
			flags_ &= ~flags;
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;

		virtual void setFn(void * fn, bool isComplex = false);
	protected:
		virtual double evaluateDRetD(double *args);
		virtual std::complex<double> evaluateDRetC(double *args);
		virtual double evaluateCRetD(std::complex<double> *args);
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args);

		// Function-pointer:
		void * fnC_;
		void * fnD_;
	};

	class MATHPRESSO_API MpOperationBinary : public MpOperation
	{
	public:
		MpOperationBinary(const Signature &s, uint32_t priority) :	MpOperation(nargs, flags),
		{
			priority_ = priority;
		}

		virtual ~MpOperationBinary()
		{
		}

		// calls generatAsmReal() and comppComplex() after setting up.
		virtual JitVar compile(JitCompiler* jc, AstNode * node) override;

		// uses calculateReal() and calculateComplex() to calculate immediate values.
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;

		bool hasFlag(uint32_t flag)
		{
			return flag & flags_;
		}

		void addFlags(uint32_t flags)
		{
			flags_ |= flags;
		}

		void removeFlags(uint32_t flags)
		{
			flags_ &= ~flags;
		}

	protected:
		// These are called by compile() and should only contain the asm-statements. vl will always
		// be in a register, vr can be in Register or in Memory.
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr);
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr);

		// Used to calculate optimization of immediates.
		virtual double calculateReal(double vl, double vr) 
		{ 
			return std::numeric_limits<double>::quiet_NaN(); 
		};
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) 
		{ 
			return std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()); 
		};
		Signature signature_;
	};
}


#endif //_MP_OPERATION_P_H