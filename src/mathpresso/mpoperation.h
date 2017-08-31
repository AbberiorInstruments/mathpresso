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

namespace mathpresso
{

	// Forward Declarations:
	struct JitCompiler;
	struct AstNode;
	struct AstOptimizer;
	struct JitVar;
	struct Context;

	typedef JitVar(*mpAsmFunc)(JitCompiler*, JitVar*);

	struct MATHPRESSO_API Signature
	{
		enum class type
		{
			real = 0,
			complex = 1
		};

		struct param
		{
			type type_;
			std::string name_;
		};

		Signature(type retType, std::vector<param> params);
		Signature(size_t nargs, type paramType = type::real);

		bool areParams(type _type) const;
		std::string to_string();

		type return_type_;
		std::vector<param> parameters_;
	private :
		std::string typeToString(type _type);
		void init(type retType, std::vector<param> params);
	};

	class MATHPRESSO_API MpOperation
	{
	public:
		enum Flags
		{
			None = 0,
			FlagHasState = 0x00000002,
			RighttoLeft = 0x00000008,
			IsAssignment = 0x0001000,
		};

		MpOperation(const Signature &s, uint32_t flags, uint32_t precedence = 0) :
			signature_(s),
			flags_(flags),
			precedence_(precedence)
		{
		}

		virtual ~MpOperation()
		{
		}

		// Add ASM code to compiler stack 
		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode>  node) const = 0;
		// Optimize AST 
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const = 0;

		size_t nargs() const
		{
			return signature_.parameters_.size();
		}
		bool isRightToLeft() const
		{
			return (flags_ & RighttoLeft) != 0;
		}
		uint32_t precedence() const
		{
			return precedence_;
		}
		uint32_t flags() const
		{
			return flags_;
		}
		Signature signature() const
		{
			return signature_;
		}
	protected:
		Signature signature_;
		uint32_t flags_;
		uint32_t precedence_;
	};

	template<typename RET, typename PARAM>
	class MATHPRESSO_API MpOperationFunc : public MpOperation
	{
	public:
		// Con-/Destructor
		MpOperationFunc(const Signature & signature, void * fnPtr, uint32_t priority = 0) : MpOperation(signature, priority),
			fnPtr_(fnPtr)
		{
		}

		// Con-/Destructor
		MpOperationFunc(uint32_t flags, size_t numargs, void * fnPtr, uint32_t priority = 0);

		virtual ~MpOperationFunc()
		{
		}

		bool hasFlag(uint32_t flag) const
		{
			return flag & flags_;
		}

		virtual JitVar compile(JitCompiler *jc, std::shared_ptr<AstNode> node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;

	protected:
		virtual RET evaluate(PARAM * args) const;

		// Function-pointer:
		void * fnPtr_;
	};

	template<typename T>
	class MATHPRESSO_API MpOperationBinary : public MpOperation
	{
	public:
		enum Flags
		{
			NopIfLZero = 0x10000000,
			NopIfRZero = 0x20000000,
			NopIfLOne = 0x40000000,
			NopIfROne = 0x80000000,
			NopIfZero = NopIfLZero | NopIfRZero,
			NopIfOne = NopIfLOne | NopIfROne,
			IsCommutativ = 0x00000010
		};

		MpOperationBinary(const Signature &signature, uint32_t flags, uint32_t priority) :
			MpOperation(signature, flags, priority)
		{
		}

		virtual ~MpOperationBinary()
		{
		}

		// calls generatAsm() after setting up.
		virtual JitVar compile(JitCompiler* jc, std::shared_ptr<AstNode> node) const override;

		// uses calculate() to calculate immediate values.
		virtual uint32_t optimize(AstOptimizer *opt, std::shared_ptr<AstNode> node) const override;

		bool hasFlag(uint32_t flag) const
		{
			return flag & flags_;
		}
	protected:
		// These are called by compile() and should only contain the asm-statements. vl will always
		// be in a register, vr can be in Register or in Memory.
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const;

		// Used to calculate optimization of immediates.
		virtual T calculate(T vl, T vr)  const;
	};

}


#endif //_MP_OPERATION_P_H