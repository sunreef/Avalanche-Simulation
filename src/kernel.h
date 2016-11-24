#pragma once

#define _USE_MATH_DEFINES 

#include <iostream>
#include <math.h>

#define INV_PI 1.0/M_PI

class Kernel
{
public:
	Kernel(float smoothingLength = 0.1f): m_invCubeSmoothingLength(1.0f / std::pow(smoothingLength,3)) {}
	~Kernel() {}

	virtual float evaluate(float q) const { return 0.0f; }
	virtual float evaluateDeriv(float q) const { return 0.0f; }

protected:
	float m_invCubeSmoothingLength;

	void printError() const {
		std::cout << "Invalid argument for kernel function" << std::endl;
	}
};



class CubicKernel : public Kernel {
public:
	CubicKernel(float smoothingLength = 0.1f): Kernel(smoothingLength) {}
	~CubicKernel() {}

	float evaluate(float q) const override {
		if (q < 0) {
			printError();
			return 0.0f;
		}
		if (q > 2) {
			return 0.0f;
		}

		float result = 0.25f * INV_PI * m_invCubeSmoothingLength;
		if (q < 1) {
			result *= 4.0f + q*q *(3.0f * q - 6.0f);
		}
		else {
			result *= std::pow(2.0f - q, 3);
		}
		return result;
	}

	float evaluateDeriv(float q) const override {
		if (q < 0) {
			printError();
			return 0.0f;
		}
		if (q > 2) {
			return 0.0f;
		}

		float result = 0.25f * INV_PI * m_invCubeSmoothingLength;
		if (q < 1) {
			result *= q*(9.0f * q - 12.0f);
		}
		else {
			result *= - 3.0f * std::pow(2.0f - q, 3);
		}
		return result;
	}
};

class ExponentialKernel : public Kernel {
public:
	ExponentialKernel(float smoothingLength = 0.1f) : Kernel(smoothingLength), sigma(1.0f / (std::pow(M_PI, 1.5)) * m_invCubeSmoothingLength) {}
	~ExponentialKernel() {}

	float evaluate(float q) const override {		
		return sigma * m_invCubeSmoothingLength * std::exp(- q*q);
	}

	float evaluateDeriv(float q) const override {
		return sigma * m_invCubeSmoothingLength * (-2 * q) * std::exp(-q*q);
	}

private:
	float sigma;
};