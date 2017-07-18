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

		void removeFlags(uint32_t flags) {
			flags_ &= ~flags;
		}

	protected:
		std::string description_ = "No description given.";
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
				addFlags(OpHasNoReal);
			}
			if (!fnC)
			{
				addFlags(OpHasNoComplex);
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
		virtual std::complex<double> evaluateDRetC(double *args);
		virtual double evaluateCRetD(std::complex<double> *args);
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args);

		// Should be overridden, if there is a special opportunity for optimization.
		// Will be called by optimize as the last operation.
		virtual uint32_t optimizeSpecial(AstOptimizer *opt, AstNode *node);

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
	
	// isfinite
	class MpOperationIsFinite :public MpOperationFuncAsm
	{
	public:
		MpOperationIsFinite() :
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal | OpHasNoComplex);
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	private:
		virtual double evaluateDRetD(double *args) override ;
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;
		
	};

	// isinf
	class MpOperationIsInfinite :public MpOperationFuncAsm
	{
	public:
		MpOperationIsInfinite() :
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal | OpHasNoComplex);
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	private:
		virtual double evaluateDRetD(double *args) override;
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;

	};

	// isnan
	class MpOperationIsNan :public MpOperationFuncAsm
	{
	public:
		MpOperationIsNan() :
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal | OpHasNoComplex);
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	private:
		virtual double evaluateDRetD(double *args) override;
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;

	};

	// real
	class MpOperationGetReal :public MpOperationFuncAsm
	{
	public:
		MpOperationGetReal() :
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoComplex);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	private:		
		virtual double evaluateCRetD(std::complex<double> *args) override;
	};

	// imag
	class MpOperationGetImag :public MpOperationFuncAsm
	{
	public:
		MpOperationGetImag() :
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoComplex);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	private:
		virtual double evaluateCRetD(std::complex<double> *args) override;
	};
	
	// Square root
	class MpOperationSqrt : public MpOperationFuncAsm
	{
	public:
		MpOperationSqrt() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		virtual double evaluateDRetD(double *args) override;

	};

	// Square root, complex result
	class MpOperationSqrtC : public MpOperationFunc
	{
	public:
		MpOperationSqrtC();
	private:
	};

	// Negation
	class MpOperationNeg : public MpOperationFuncAsm
	{
	public:
		MpOperationNeg() : MpOperationFuncAsm(1, OpFlagHasAsm, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal | OpHasNoComplex);
			priority_ = 3;
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	protected:
		virtual double evaluateDRetD(double *args) override;
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;
		virtual uint32_t optimizeSpecial(AstOptimizer *opt, AstNode *node) override;
	};

	// Not
	class MpOperationNot : public MpOperationFuncAsm
	{
	public:
		MpOperationNot() : MpOperationFuncAsm(1, OpFlagHasAsm, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoReal | OpHasNoComplex);
			priority_ = 3;
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	protected:
		virtual double evaluateDRetD(double *args) override;
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;
	};

	// conjugate
	class MpOperationConjug :public MpOperationFuncAsm
	{
	public:
		MpOperationConjug() :
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoComplex);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	private:
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;
		virtual uint32_t optimizeSpecial(AstOptimizer *opt, AstNode *node) override;
	};

	// Reciprocate
	class MpOperationRecip :public MpOperationFuncAsm
	{
	public:
		MpOperationRecip() :
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoComplex | OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	private:
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;
		virtual double evaluateDRetD(double * args) override;
	};

	// trigonometric functions (sin, cos, tan, etc)
	class MpOperationTrigonometrie : public MpOperationFunc
	{
	public:
		enum trigonometrieFunc {
			sin, cos, tan, asin, acos, atan, sinh, cosh, tanh
		};

		MpOperationTrigonometrie(uint32_t type) : 
			MpOperationFunc(1, OpFlagHasAsm, nullptr, nullptr),
			type_(type)
		{
			removeFlags(OpHasNoReal | OpHasNoComplex);
		}
		virtual ~MpOperationTrigonometrie() {}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	protected:
		virtual double evaluateDRetD(double *args) override;
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;
	private:
		uint32_t type_;
	};
	
	// Sign bit
	class MpOperationSignBit : public MpOperationFuncAsm
	{
	public:
		MpOperationSignBit() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		virtual double evaluateDRetD(double *args) override;
	};

	// copy sign
	class MpOperationCopySign : public MpOperationFuncAsm
	{
	public:
		MpOperationCopySign() : MpOperationFuncAsm(2, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		virtual double evaluateDRetD(double *args) override;
	};

	// Average
	class MpOperationAvg : public MpOperationFuncAsm
	{
	public:
		MpOperationAvg() : MpOperationFuncAsm(2, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr) 
		{
			removeFlags(OpHasNoComplex | OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	private:
		virtual double evaluateDRetD(double *args) override;
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;
	};

	// Absolute
	class MpOperationAbs : public MpOperationFuncAsm
	{
	public:
		MpOperationAbs() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoComplex | OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	private:
		virtual double evaluateDRetD(double *args) override;
		virtual double evaluateCRetD(std::complex<double> *args) override;
	};

	// round
	class MpOperationRound : public MpOperationFuncAsm
	{
	public: 
		MpOperationRound() : MpOperationFuncAsm(1, OpFlagNone, nullptr, nullptr, nullptr, nullptr) 	
		{
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		std::string description_ = "Rounds *.5 towards infinity.";
		virtual double evaluateDRetD(double *args) override;
	};

	// roundeven
	class MpOperationRoundEven : public MpOperationFuncAsm
	{
	public:
		MpOperationRoundEven() : MpOperationFuncAsm(1, OpFlagNone, nullptr, nullptr, nullptr, nullptr) {
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		std::string description_ = "Rounds *.5 towards nearest even integer.";
		virtual double evaluateDRetD(double *args) override;
	};

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
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr);
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr);

		// Used to calculate optimization of immediates.
		virtual double calculateReal(double vl, double vr) { return NAN; };
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) { return std::complex<double>(NAN, NAN); };
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

	// Minimum
	class MpOperationMin : public MpOperationBinary
	{
	public:
		MpOperationMin() :
			MpOperationBinary(2, MpOperationFlags::OpFlagHasAsm, 0) {
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
	};

	// Maximum
	class MpOperationMax : public MpOperationBinary
	{
	public:
		MpOperationMax() :
			MpOperationBinary(2, MpOperationFlags::OpFlagHasAsm, 0) {
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
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
		virtual double calculateReal(double vl, double vr) override;
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
		virtual double calculateReal(double vl, double vr) override;
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
		virtual double calculateReal(double vl, double vr) override;
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
		virtual double calculateReal(double vl, double vr) override;
	};

	// Modulo
	class MpOperationModulo : public MpOperationBinary
	{
	public:
		MpOperationModulo() :
			MpOperationBinary(2, MpOperationFlags::OpHasNoComplex | MpOperationFlags::OpFlagHasAsm, 5) {
		}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
	};

	// Ternary Operation
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

	// Assignment
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