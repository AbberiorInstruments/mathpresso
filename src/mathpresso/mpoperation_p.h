#pragma once
#include <mathpresso/mpoperation.h>

namespace mathpresso
{
	// isfinite
	class MpOperationIsFinite : public MpOperationFunc
	{
	public:
		MpOperationIsFinite();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// isinf
	class MpOperationIsInfinite : public MpOperationFunc
	{
	public:
		MpOperationIsInfinite();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// isnan
	class MpOperationIsNan : public MpOperationFunc
	{
	public:
		MpOperationIsNan();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// real
	class MpOperationGetReal : public MpOperationFunc
	{
	public:
		MpOperationGetReal();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// imag
	class MpOperationGetImag : public MpOperationFunc
	{
	public:
		MpOperationGetImag();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Square root
	class MpOperationSqrt : public MpOperationFunc
	{
	public:
		MpOperationSqrt();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Negation
	class MpOperationNeg : public MpOperationFunc
	{
	public:
		MpOperationNeg();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	protected:
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;
	};

	// Not
	class MpOperationNot : public MpOperationFunc
	{
	public:
		MpOperationNot();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// conjugate
	class MpOperationConjug :public MpOperationFunc
	{
	public:
		MpOperationConjug();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	private:
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;
	};

	// Reciprocate
	class MpOperationRecip :public MpOperationFunc
	{
	public:
		MpOperationRecip();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Sign bit
	class MpOperationSignBit : public MpOperationFunc
	{
	public:
		MpOperationSignBit();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// copy sign
	class MpOperationCopySign : public MpOperationFunc
	{
	public:
		MpOperationCopySign();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Average
	class MpOperationAvg : public MpOperationFunc
	{
	public:
		MpOperationAvg();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Absolute
	class MpOperationAbs : public MpOperationFunc
	{
	public:
		MpOperationAbs();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// round
	class MpOperationRound : public MpOperationFunc
	{
	public:
		MpOperationRound();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	protected:
		std::string description_ = "Rounds *.5 towards infinity.";
	};

	// roundeven
	class MpOperationRoundEven : public MpOperationFunc
	{
	public:
		MpOperationRoundEven();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	protected:
		std::string description_ = "Rounds *.5 towards nearest even integer.";
	};

	// trunc
	class MpOperationTrunc : public MpOperationFunc
	{
	public:
		MpOperationTrunc();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// frac
	class MpOperationFrac : public MpOperationFunc
	{
	public:
		MpOperationFrac();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// floor
	class MpOperationFloor : public MpOperationFunc
	{
	public:
		MpOperationFloor();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// ceil
	class MpOperationcCeil : public MpOperationFunc
	{
	public:
		MpOperationcCeil();
		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
	};

	// Addition
	template<typename T>
	class MpOperationAdd : public MpOperationBinarytemp<T>
	{
	public:
		MpOperationAdd();

	protected:

		virtual JitVar generateAsm(JitCompiler * jc, JitVar vl, JitVar vr) const override;

		virtual T calculate(T vl, T vr) const override;
	};

	// Subtraction
	class MpOperationSub : public MpOperationBinary
	{
	public:
		MpOperationSub() :
			MpOperationBinary(Signature(2, Signature::type::both, MpOperationFlags::OpFlagNopIfRZero  | MpOperationFlags::OpFlagIsOperator), 6)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) const override;
	};

	// Multiplication
	class MpOperationMul : public MpOperationBinary
	{
	public:
		MpOperationMul() :
			MpOperationBinary(Signature(2, Signature::type::both, MpOperationFlags::OpFlagNopIfZero | MpOperationFlags::OpIsCommutativ  | MpOperationFlags::OpFlagIsOperator), 5)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) const override;
	};

	// Division
	class MpOperationDiv : public MpOperationBinary
	{
	public:
		MpOperationDiv() :
			MpOperationBinary(Signature(2, Signature::type::both, MpOperationFlags::OpFlagNopIfLOne  | MpOperationFlags::OpFlagIsOperator), 5)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) const override;
	};

	// Minimum
	class MpOperationMin : public MpOperationBinary
	{
	public:
		MpOperationMin() :
			MpOperationBinary(Signature(2), 0)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
	};

	// Maximum
	class MpOperationMax : public MpOperationBinary
	{
	public:
		MpOperationMax() :
			MpOperationBinary(Signature(2), 0)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
	};

	// Equality
	class MpOperationEq : public MpOperationBinary
	{
	public:
		MpOperationEq() :
			MpOperationBinary(Signature(2, Signature::type::both, MpOperationFlags::OpIsCommutativ  | MpOperationFlags::OpFlagIsOperator), 9)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) const override;
	};

	// Inequality
	class MpOperationNe : public MpOperationBinary
	{
	public:
		MpOperationNe() :
			MpOperationBinary(Signature(2, Signature::type::both, MpOperationFlags::OpIsCommutativ  | MpOperationFlags::OpFlagIsOperator), 9)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual JitVar generateAsmComplex(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
		virtual std::complex<double> calculateComplex(std::complex<double> vl, std::complex<double> vr) const override;
	};

	// Lesser than
	class MpOperationLt : public MpOperationBinary
	{
	public:
		MpOperationLt() :
			MpOperationBinary(Signature(2, Signature::type::real, MpOperationFlags::OpHasNoComplex  | MpOperationFlags::OpFlagIsOperator), 8)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
	};

	// Lesser Equal
	class MpOperationLe : public MpOperationBinary
	{
	public:
		MpOperationLe() :
			MpOperationBinary(Signature(2, Signature::type::real, MpOperationFlags::OpHasNoComplex  | MpOperationFlags::OpFlagIsOperator), 8)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
	};

	// Greater than
	class MpOperationGt : public MpOperationBinary
	{
	public:
		MpOperationGt() :
			MpOperationBinary(Signature(2, Signature::type::real, MpOperationFlags::OpHasNoComplex  | MpOperationFlags::OpFlagIsOperator), 8)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
	};

	// Greater equal
	class MpOperationGe : public MpOperationBinary
	{
	public:
		MpOperationGe() :
			MpOperationBinary(Signature(2, Signature::type::real, MpOperationFlags::OpHasNoComplex  | MpOperationFlags::OpFlagIsOperator), 8)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
	};

	// Modulo
	class MpOperationModulo : public MpOperationBinary
	{
	public:
		MpOperationModulo() :
			MpOperationBinary(Signature(2, Signature::type::real, MpOperationFlags::OpHasNoComplex  | MpOperationFlags::OpFlagIsOperator), 5)
		{
		}

	protected:
		virtual JitVar generateAsmReal(JitCompiler * jc, JitVar vl, JitVar vr) const override;
		virtual double calculateReal(double vl, double vr) const override;
	};

	// Ternary Operation
	class MpOperationTernary : public MpOperation
	{
	public:
		MpOperationTernary() : MpOperation(Signature(3, Signature::type::both, MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpFlagIsOperator), 15)
		{
		}

		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;
	protected:
		bool isColon_;
	};

	// Assignment
	class MpOperationAssignment : public MpOperation
	{
	public:
		MpOperationAssignment() : MpOperation(Signature(1, Signature::type::both, MpOperationFlags::OpIsRighttoLeft | MpOperationFlags::OpIsAssgignment), 15)
		{
		}


		virtual JitVar compile(JitCompiler *jc, AstNode *node) const override;
		virtual uint32_t optimize(AstOptimizer *opt, AstNode *node) const override;
	};
}