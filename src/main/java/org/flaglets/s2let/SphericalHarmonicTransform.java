
package org.flaglets.s2let;

import org.bridj.*;
import java.util.Random;
import static org.flaglets.s2let.bindings.S2letLibrary.*;

/**
 *
 * @author bl
 */
public class SphericalHarmonicTransform {
    private Pointer<ComplexDouble> coefficients;
    private int bandlimit;
    private boolean reality;
    public Pointer<ComplexDouble> getCoefficients() {
        return coefficients;
    }
    public int getBandlimit() {
        return this.bandlimit;
    }
    public boolean getReality() {
        return this.reality;
    }
    public void detectReality() {
        Pointer<Double> pd = coefficients.as(Double.class);
    	for (int el = 0; el < bandlimit; el++) {
            int em = 0;
            int i = el*el + el + em;
            double val = pd.getDoubleAtOffset(i * 2 * 8 + 8);
            if( val != 0.0 ){
                this.reality = false;
                return;
            }
            for (em = 1; em <= el; em++) {
                i = el*el + el + em;
                double val1 = pd.getDoubleAtOffset(i * 2 * 8);
                double val2 = pd.getDoubleAtOffset(i * 2 * 8 + 8);
                int iop = el*el + el - em;
                double val3 = pd.getDoubleAtOffset(iop * 2 * 8);
                double val4 = pd.getDoubleAtOffset(iop * 2 * 8 + 8);
                if( val3 != Math.pow(-1.0,em) * val1 || val4 != - Math.pow(-1.0,em+1) * val2 ){
                    this.reality = false;
                    return;
                }
            }
            this.reality = true;
        }
    }
    public SphericalHarmonicTransform(int bandlimit) {
        this.coefficients = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        this.bandlimit = bandlimit;
    }
    public void randomCoefficiens() {
	Random randomGenerator = new Random();
	Pointer<Double> pd = coefficients.as(Double.class);
    	for (int i = 0; i < bandlimit*bandlimit; i++) {
            pd.setDoubleAtOffset(i * 2 * 8, randomGenerator.nextDouble());
            pd.setDoubleAtOffset(i * 2 * 8 + 8, randomGenerator.nextDouble());
        }
        this.detectReality();
    }
    public void randomCoefficiensReal() {
	Random randomGenerator = new Random();
	Pointer<Double> pd = coefficients.as(Double.class);
    	for (int el = 0; el < bandlimit; el++) {
            int em = 0;
            int i = el*el + el + em;
            pd.setDoubleAtOffset(i * 2 * 8, randomGenerator.nextDouble());
            pd.setDoubleAtOffset(i * 2 * 8 + 8, 0.0);
            for (em = 1; em <= el; em++) {
                i = el*el + el + em;
                pd.setDoubleAtOffset(i * 2 * 8, randomGenerator.nextDouble());
                pd.setDoubleAtOffset(i * 2 * 8 + 8, randomGenerator.nextDouble());
                int iop = el*el + el - em;
                pd.setDoubleAtOffset(iop * 2 * 8, Math.pow(-1.0,em) * pd.getDoubleAtOffset(i * 2 * 8));
                pd.setDoubleAtOffset(iop * 2 * 8 + 8, Math.pow(-1.0,em+1) * pd.getDoubleAtOffset(i * 2 * 8 + 8));
            }
        }
        this.detectReality();
    }
    public double maxAbsoluteDifferenceWith(Pointer<ComplexDouble> f2) {
        Pointer<ComplexDouble> f1 = coefficients;
        double val = 0.0;
        Pointer<Double> v1 = f1.as(Double.class);
        Pointer<Double> v2 = f2.as(Double.class);
        int len = Math.min( (int)f1.getValidElements(), (int)f2.getValidElements());
        for (int i = 0; i < len; i++) {
            val = Math.max( val, Math.max( v1.getDoubleAtOffset(i * 2 * 8) - v2.getDoubleAtOffset(i * 2 * 8), 
                    v1.getDoubleAtOffset(i * 2 * 8 + 8) - v2.getDoubleAtOffset(i * 2 * 8 + 8) ) );
        }
        return val;
    }
    public double maxAbsoluteDifferenceWith(SphericalHarmonicTransform f2) {
        return maxAbsoluteDifferenceWith(f2.getCoefficients());
    }
    public SphericalHarmonicTransform(Pointer<ComplexDouble> flm) {
        this.coefficients = flm;
        this.bandlimit = (int)Math.sqrt(flm.getValidElements());
        this.detectReality();
    }
    public SphericalHarmonicTransform(MWCmplxMap f) {
        this.bandlimit = f.getResolution();
        this.reality = false;
        this.coefficients = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        s2let_mw_map2alm(this.coefficients, f.getMap(), f.getResolution());
    }
    public SphericalHarmonicTransform(MWRealMap f) {
        this.bandlimit = f.getResolution();
        this.reality = true;
        this.coefficients = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        s2let_mw_map2alm_real(this.coefficients, f.getMap(), f.getResolution());
    }
    public SphericalHarmonicTransform(HealpixRealMap f, int bandlimit) {
        this.reality = true;
        this.coefficients = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        s2let_hpx_map2alm_real(this.coefficients, f.getMap(), f.getResolution(), bandlimit);
    }
    public SphericalHarmonicTransform(AxisymmetricWaveletTransform wavelets) {
        this.reality = wavelets.getReality();
        // Calculate wavelet decomposition
    }
    
}
