#pragma once
#include <mathpresso/mpoperation.h>

namespace mathpresso
{
	// isfinite
	class MpOperationIsFinite : public MpOperationFunc
	{
	public:
		MpOperationIsFinite();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// isinf
	class MpOperationIsInfinite : public MpOperationFunc
	{
	public:
		MpOperationIsInfinite();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// isnan
	class MpOperationIsNan : public MpOperationFunc
	{
	public:
		MpOperationIsNan();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// real
	class MpOperationGetReal : public MpOperationFunc
	{
	public:
		MpOperationGetReal();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// imag
	class MpOperationGetImag : public MpOperationFunc
	{
	public:
		MpOperationGetImag();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// Square root
	class MpOperationSqrt : public MpOperationFunc
	{
	public:
		MpOperationSqrt();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// Negation
	class MpOperationNeg : public MpOperationFunc
	{
	public:
		MpOperationNeg();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	protected:
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;
	};

	// Not
	class MpOperationNot : public MpOperationFunc
	{
	public:
		MpOperationNot();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// conjugate
	class MpOperationConjug :public MpOperationFunc
	{
	public:
		MpOperationConjug();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	private:
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;
	};

	// Reciprocate
	class MpOperationRecip :public MpOperationFunc
	{
	public:
		MpOperationRecip();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// trigonometric functions (sin, cos, tan, etc)
	class MpOperationTrigonometrie : public MpOperationFunc
	{
	public:
		enum class Type {
			sin, cos, tan, asin, acos, atan, sinh, cosh, tanh
		};
		MpOperationTrigonometrie(Type type);
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	protected:
		virtual double evaluateDRetD(double *args) override;
		virtual std::complex<double> evaluateCRetC(std::complex<double> *args) override;
	private:
		Type type_;
	};

	// Sign bit
	class MpOperationSignBit : public MpOperationFunc
	{
	public:
		MpOperationSignBit();		
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// copy sign
	class MpOperationCopySign : public MpOperationFunc
	{
	public:
		MpOperationCopySign();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// Average
	class MpOperationAvg : public MpOperationFunc
	{
	public:
		MpOperationAvg();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// Absolute
	class MpOperationAbs : public MpOperationFunc
	{
	public:
		MpOperationAbs();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// round
	class MpOperationRound : public MpOperationFunc
	{
	public:
		MpOperationRound();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	protected:
		std::string description_ = "Rounds *.5 towards infinity.";
	};

	// roundeven
	class MpOperationRoundEven : public MpOperationFunc
	{
	public:
		MpOperationRoundEven();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	protected:
		std::string description_ = "Rounds *.5 towards nearest even integer.";
	};

	// trunc
	class MpOperationTrunc : public MpOperationFunc
	{
	public:
		MpOperationTrunc();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// frac
	class MpOperationFrac : public MpOperationFunc
	{
	public:
		MpOperationFrac();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// floor
	class MpOperationFloor : public MpOperationFunc
	{
	public:
		MpOperationFloor();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
	};

	// ceil
	class MpOperationcCeil : public MpOperationFunc
	{
	public:
		MpOperationcCeil();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
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
		MpOperationTernary(bool iscolon) : MpOperation(2, MpOperationFlags::OpIsRighttoLeft), isColon_(iscolon)
		{
			priority_ = 15;
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;
	protected:
		bool isColon_;
	};

	// Assignment
	class MpOperationAssignment : public MpOperation
	{
	public:
		MpOperationAssignment() : MpOperation(2, MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpIsAssgignment)
		{
			priority_ = 15;
		}


		virtual JitVar compile(JitCompiler *jc, AstNode *node) override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) override;
	};
}