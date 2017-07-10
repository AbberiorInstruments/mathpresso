// [MathPresso]
// Mathematical Expression Parser and JIT Compiler.
//
// [License]
// Zlib - See LICENSE.md file in the package.

// [Guard]
#ifndef _MP_OPERATION_P_H
#define _MP_OPERATION_P_H

#include <complex>


namespace mathpresso {

	// Forward Declarations:
	struct JitCompiler;
	struct AstNode;
	struct AstOptimizer;
	struct JitVar;


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
		
		// Set for optimization of binary operations
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
			flags_(flags),
			priority_(0)
		{
		}

		virtual ~MpOperation() 
		{
		}

		// Add ASM code to compiler stack 
		virtual JitVar compile(JitCompiler *jc, AstNode * node)  = 0;
		// Optimize AST 
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) = 0;

		uint32_t flags() 
		{ 
			return flags_; 
		}
		uint32_t nargs() 
		{ 
			return nargs_; 
		}

		bool hasFlag(uint32_t flag) 
		{
			return flag & flags_;
		}

		void addFlags(uint32_t flags)
		{
			flags_ |= flags;
		}

	protected:
		uint32_t nargs_;
		uint32_t flags_;
		uint32_t priority_;
	};

	class MpOperationFunc : public MpOperation
	{
	public:
		// Con-/Destructor
		MpOperationFunc(uint32_t nargs, uint32_t flags, void * fnD, void * fnC) : MpOperation(nargs, flags),
			fnD_(fnD),
			fnC_(fnC)
		{
			if (!fnD) 
			{
				flags_ |= OpHasNoReal;
			}
			if (!fnC)
			{
				flags_ |= OpHasNoComplex;
			}
		}

		virtual ~MpOperationFunc()
		{
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;


		virtual void setFn(void * fn, bool isComplex = false);
	protected:
		virtual double evaluateDRetD(double *args);
		
		double fnDouble(AstOptimizer * opt, AstNode *node);
		std::complex<double> fnComplex(AstOptimizer * opt, AstNode *node);

		// Function-pointer:
		void * fnC_;
		void * fnD_;
	};

	class MpOperationFuncAsm : public MpOperationFunc
	{
	public:
		MpOperationFuncAsm(uint32_t nargs, uint32_t flags, void * fnD, void * fnC, mpAsmFunc asmC, mpAsmFunc asmD) :
			MpOperationFunc(nargs, flags | MpOperationFlags::OpFlagHasAsm, fnD, fnC),
			asmC_(asmC),
			asmD_(asmD)
		{
		}

		virtual ~MpOperationFuncAsm()
		{
		}
		
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

		virtual void setFnAsm(mpAsmFunc fn, bool isComplex = false);

	protected:
		mpAsmFunc asmC_;
		mpAsmFunc asmD_;
	};


	//class MpOprationUnary :public MpOperationFuncAsm
	//{
	//public:
	//	//MpOprationUnary() :
	//	//	MpOperationFuncAsm(1, 0, (void*)calculateReal, (void)calculateComplex, generateAsmReal, generateAsmComplex)
	//	//{
	//	//}

	//	virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;
	//private:
	//	JitVar generateAsmComplex(JitCompiler *jc, AstNode * node);
	//	 JitVar generateAsmReal(JitCompiler *jc, AstNode * node);
	//	double calculateReal(double* args);
	//	std::complex<double> calculateComplex(std::complex<double>* args);
	//};

	// ============================================================================
	// Binary operations
	// ============================================================================
	class MpOperationBinary : public MpOperation
	{
	public:
		MpOperationBinary(uint32_t nargs, uint32_t  flags, uint32_t priority) :
			MpOperation(nargs, flags)
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

	protected:
		// These are called by compile() and should only contain the asm-statements.
		// vl will always be in a register, vr can be in Register or in Memory.
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) = 0;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) = 0;

		// Used to calculate optimization of immediates.
		virtual double calculateReal(double vl, double vr) = 0;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) = 0;
	};
		
	// Addition
	class MpOperationAdd : public MpOperationBinary
	{
	public:
		MpOperationAdd() :
			MpOperationBinary(2, MpOperationFlags::OpFlagNopIfZero | MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagHasAsm, 6)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Subtraction
	class MpOperationSub : public MpOperationBinary
	{
	public:
		MpOperationSub() :
			MpOperationBinary(2, MpOperationFlags::OpFlagNopIfRZero | MpOperationFlags::OpFlagHasAsm, 6)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Multiplication
	class MpOperationMul : public MpOperationBinary
	{
	public:
		MpOperationMul() :
			MpOperationBinary(2, MpOperationFlags::OpFlagNopIfZero | MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagHasAsm, 5)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Division
	class MpOperationDiv : public MpOperationBinary
	{
	public:
		MpOperationDiv() :
			MpOperationBinary(2, MpOperationFlags::OpFlagNopIfLOne | MpOperationFlags::OpFlagHasAsm, 5)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Equality
	class MpOperationEq : public MpOperationBinary
	{
	public:
		MpOperationEq() :
			MpOperationBinary(2, MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagHasAsm, 9)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Inequality
	class MpOperationNe : public MpOperationBinary
	{
	public:
		MpOperationNe() :
			MpOperationBinary(2, MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagHasAsm, 9)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Lesser than
	class MpOperationLt : public MpOperationBinary
	{
	public:
		MpOperationLt() :
			MpOperationBinary(2, MpOperationFlags::OpHasNoComplex | MpOperationFlags::OpFlagHasAsm, 8)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Lesser Equal
	class MpOperationLe : public MpOperationBinary
	{
	public:
		MpOperationLe() :
			MpOperationBinary(2, MpOperationFlags::OpHasNoComplex | MpOperationFlags::OpFlagHasAsm, 8)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Greater than
	class MpOperationGt : public MpOperationBinary
	{
	public:
		MpOperationGt() :
			MpOperationBinary(2, MpOperationFlags::OpHasNoComplex | MpOperationFlags::OpFlagHasAsm, 8)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	// Greater equal
	class MpOperationGe : public MpOperationBinary
	{
	public:
		MpOperationGe() :
			MpOperationBinary(2, MpOperationFlags::OpHasNoComplex | MpOperationFlags::OpFlagHasAsm, 8)
		{
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) override;
	};

	class MpOperationTernary : public MpOperation
	{
	public:
		MpOperationTernary() : MpOperation(3, MpOperationFlags::OpIsRighttoLeft)
		{
			priority_ = 15;
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;

	};

	class MpOperationAssignment : public MpOperation
	{
	public:
		MpOperationAssignment() : MpOperation(2, MpOperationFlags::OpIsRighttoLeft) 
		{
			priority_ = 15; 
		}


		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;
	};

}


#endif //_MP_OPERATION_P_H