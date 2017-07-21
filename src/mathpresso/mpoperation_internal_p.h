#pragma once
#include "mpoperation_p.h"

namespace mathpresso
{
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
		virtual double evaluateDRetD(double *args) override;
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
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, nullptr, nullptr, nullptr, nullptr)
		{
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
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, nullptr, nullptr, nullptr, nullptr)
		{
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
		MpOperationNot() : MpOperationFuncAsm(1, OpFlagHasAsm, nullptr, nullptr, nullptr, nullptr)
		{
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
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
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
			MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
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
		MpOperationSignBit() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
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
		MpOperationCopySign() : MpOperationFuncAsm(2, MpOperationFlags::OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
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
		MpOperationAbs() : MpOperationFuncAsm(1, MpOperationFlags::OpFlagCReturnsD, nullptr, nullptr, nullptr, nullptr)
		{
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
		MpOperationRoundEven() : MpOperationFuncAsm(1, OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		std::string description_ = "Rounds *.5 towards nearest even integer.";
		virtual double evaluateDRetD(double *args) override;
	};

	// trunc
	class MpOperationTrunc : public MpOperationFuncAsm
	{
	public:
		MpOperationTrunc() : MpOperationFuncAsm(1, OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		virtual double evaluateDRetD(double *args) override;
	};

	// frac
	class MpOperationFrac : public MpOperationFuncAsm
	{
	public:
		MpOperationFrac() : MpOperationFuncAsm(1, OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		virtual double evaluateDRetD(double *args) override;
	};

	// floor
	class MpOperationFloor : public MpOperationFuncAsm
	{
	public:
		MpOperationFloor() : MpOperationFuncAsm(1, OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		virtual double evaluateDRetD(double *args) override;
	};

	// ceil
	class MpOperationcCeil : public MpOperationFuncAsm
	{
	public:
		MpOperationcCeil() : MpOperationFuncAsm(1, OpFlagNone, nullptr, nullptr, nullptr, nullptr)
		{
			removeFlags(OpHasNoReal);
		}
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;

	protected:
		virtual double evaluateDRetD(double *args) override;
	};

	// Log
	class MpOperationLog : public MpOperationFunc
	{
	public:
		MpOperationLog();
	};

	// Log2
	class MpOperationLog2 : public MpOperationFunc
	{
	public:
		MpOperationLog2();
	};

	// Log10
	class MpOperationLog10 : public MpOperationFunc
	{
	public:
		MpOperationLog10();
	};

	// exp
	class MpOperationExp : public MpOperationFunc
	{
	public:
		MpOperationExp();
	private:
	};

	// pow
	class MpOperationPow : public MpOperationFunc
	{
	public:
		MpOperationPow();
	};

	// atan2
	class MpOperationAtan2 : public MpOperationFunc
	{
	public:
		MpOperationAtan2();
	};

	// hypot
	class MpOperationHypot : public MpOperationFunc
	{
	public:
		MpOperationHypot();
	};

	
	// Addition
	class MpOperationAdd : public MpOperationBinary
	{
	public:
		MpOperationAdd() :
			MpOperationBinary(2, MpOperationFlags::OpFlagNopIfZero | MpOperationFlags::OpIsCommutativ | MpOperationFlags::OpFlagHasAsm, 6)
		{}

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
		{}

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
		{}

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
		{}

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
			MpOperationBinary(2, MpOperationFlags::OpFlagHasAsm, 0)
		{}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
	};

	// Maximum
	class MpOperationMax : public MpOperationBinary
	{
	public:
		MpOperationMax() :
			MpOperationBinary(2, MpOperationFlags::OpFlagHasAsm, 0)
		{}

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
		{}

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
		{}

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
		{}

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
		{}

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
		{}

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
		{}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
	};

	// Modulo
	class MpOperationModulo : public MpOperationBinary
	{
	public:
		MpOperationModulo() :
			MpOperationBinary(2, MpOperationFlags::OpHasNoComplex | MpOperationFlags::OpFlagHasAsm, 5)
		{}

	protected:
		virtual JitVar generatAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) override;
		virtual double calculateReal(double vl, double vr) override;
	};

	// Ternary Operation
	class MpOperationTernary : public MpOperation
	{
	public:
		MpOperationTernary() : MpOperation(2, MpOperationFlags::OpIsRighttoLeft)
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