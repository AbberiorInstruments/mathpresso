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

	enum MpOperationFlags
	{
		OpFlagNone = 0,
		OpFlagIsOperator = 0x00000001,
		OpFlagHasState = 0x00000002,
		OpIsRighttoLeft = 0x00000008,

		OpIsCommutativ = 0x00000010,

		// some information, that is necessary for parsing:
		OpIsAssgignment = 0x0001000,

		// Set for optimization of binary operations
		OpFlagNopIfLZero = 0x10000000,
		OpFlagNopIfRZero = 0x20000000,
		OpFlagNopIfLOne = 0x40000000,
		OpFlagNopIfROne = 0x80000000,
		OpFlagNopIfZero = OpFlagNopIfLZero | OpFlagNopIfRZero,
		OpFlagNopIfOne = OpFlagNopIfLOne | OpFlagNopIfROne
	};

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

		Signature(type retType, std::vector<param> params)
		{
			init(retType, params);
		}

		Signature(size_t nargs, type paramType = type::real)
		{
			init(paramType, { nargs, { paramType, "" } });
		}

		bool areParams(type _type) const
		{
			bool ret = true;
			for (auto p : parameters_)
			{
				ret &= (p.type_ == _type);
			}
			return ret;
		}

		std::string to_string()
		{
			std::string out("");
			out += typeToString(return_type_);
			out += " (";
			for (size_t i = 0; i < parameters_.size(); i++)
			{
				out += typeToString(parameters_[i].type_);
				if (parameters_[i].name_ != "")
					out += " " + parameters_[i].name_;
				if (i != parameters_.size() -1)
				out += ", ";
			}
			out += ")";
			return out;
		}


		type return_type_;
		std::vector<param> parameters_;
	private :
		std::string typeToString(type _type)
		{
			switch (_type)
			{
				case type::real: return "real";
				case type::complex: return "complex";
				default:
					throw std::runtime_error("unknown type.");
			}
		}

		void init(type retType, std::vector<param> params);
	};

	class MATHPRESSO_API MpOperation
	{
	public:
		// Con-/Destructor
		MpOperation(size_t nargs, uint32_t flags, uint32_t priority = 0) :
			signature_(nargs, Signature::type::real),
			nargs_(nargs),
			flags_(flags),
			priority_(priority)
		{
		}

		MpOperation(const Signature &s, uint32_t flags, uint32_t priority = 0) :
			signature_(s),
			nargs_(s.parameters_.size()),
			flags_(flags),
			priority_(priority)
		{
		}

		virtual ~MpOperation()
		{
		}

		// Add ASM code to compiler stack 
		virtual JitVar compile(JitCompiler *jc, AstNode * node) const = 0;
		// Optimize AST 
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const = 0;

		size_t nargs() const
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
		Signature signature() const
		{
			return signature_;
		}

	protected:
		Signature signature_;
		size_t nargs_;
		uint32_t flags_;
		uint32_t priority_;
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

		void addFlags(uint32_t flags)
		{
			flags_ |= flags;
		}

		void removeFlags(uint32_t flags)
		{
			flags_ &= ~flags;
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;

	protected:
		virtual RET evaluate(PARAM * args) const;

		// Function-pointer:
		void * fnPtr_;
	};

	template<typename T>
	class MATHPRESSO_API MpOperationBinary : public MpOperation
	{
	public:
		MpOperationBinary(const Signature &signature, uint32_t flags, uint32_t priority) :
			MpOperation(signature, flags, priority)
		{
		}

		virtual ~MpOperationBinary()
		{
		}

		// calls generatAsm() after setting up.
		virtual JitVar compile(JitCompiler* jc, AstNode * node) const override;

		// uses calculate() to calculate immediate values.
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;

		bool hasFlag(uint32_t flag) const
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
		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const;

		// Used to calculate optimization of immediates.
		virtual T calculate(T vl, T vr)  const;
	};

}


#endif //_MP_OPERATION_P_H