module;

#include "tw_includes.h"
#include "tw_test.h"

/// Test of FFT routines, to compare with numpy FFT routines, we have to account for difference
/// in phase convention, which leads to interchange of positive/negative frequencies.
export module fft_test;
import base;
import tw_iterator;
import driver;
import fft;
import fields;

export struct FFTTest: ComputeTool {
    FFTTest(const std::string& name,MetricSpace *ms,Task *tsk): ComputeTool(name,ms,tsk) {}
    virtual void RegisterTests() {
        REGISTER(FFTTest,ComplexFFTTest);
        REGISTER(FFTTest,HermiticityTest);
        REGISTER(FFTTest,RealFFTTest);
    }
	void ComplexFFTTest() {
        std::array<tw::Complex,4> x;
        std::array<tw::Complex,4> k;

        x = {1,0,0,0};
        k = {1,1,1,1};
        fft::ComplexFFT((tw::Float*)&x[0],(tw::Float*)&x[0]+1,4,2,1.0);
        for (auto i=0; i<4; i++) {
            ASSERT_NEAR(x[i].real(),k[i].real(),1e-4);
            ASSERT_NEAR(x[i].imag(),k[i].imag(),1e-4);
        }

        x = {1.0 + ii,2.0 - ii,2.0*ii,0.5};
        k = {3.5 + 2.0*ii,2.0 + 0.5*ii,-1.5 + 4.0*ii,-2.5*ii};
        fft::ComplexFFT((tw::Float*)&x[0],(tw::Float*)&x[0]+1,4,2,1.0);
        for (auto i=0; i<4; i++) {
            ASSERT_NEAR(x[i].real(),k[i].real(),1e-4);
            ASSERT_NEAR(x[i].imag(),k[i].imag(),1e-4);
        }
    }
    /// Transform some real data, we should get hermiticity
    void HermiticityTest() {
        std::array<tw::Complex,8> x;
        std::array<tw::Complex,8> k;
        x = {0,0,1,2,1,0,0,0};
        k = {4,2.4142*(-1.0+ii),-2.0*ii,0.4142*(1.0+ii),0,0.4142*(1.0-ii),2.0*ii,2.4142*(-1.0-ii)};
        fft::ComplexFFT((tw::Float*)&x[0],(tw::Float*)&x[0]+1,8,2,1.0);
        for (auto i=0; i<8; i++) {
            ASSERT_NEAR(x[i].real(),k[i].real(),1e-4);
            ASSERT_NEAR(x[i].imag(),k[i].imag(),1e-4);
        }        
    }
	void RealFFTTest() {
        std::array<tw::Float,8> x;
        std::array<tw::Float,8> k; // k1,kN,k2r,k2i,...

        x = {1,0,0,0};
        k = {1,1,1,0};
        fft::RealFFT(&x[0],4,1,1);
        for (auto i=0; i<4; i++) {
            ASSERT_NEAR(x[i],k[i],1e-4);
        }

        x = {0,0,1,2,1,0,0,0};
        k = {4,0,-2.4142,2.4142,0,-2,0.41421,0.41421};
        fft::RealFFT(&x[0],8,1,1);
        for (auto i=0; i<8; i++) {
            ASSERT_NEAR(x[i],k[i],1e-4);
        }
    }
};
