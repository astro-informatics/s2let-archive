
package org.s2let;

import org.bridj.*;
import java.util.Random;

/**
 *
 * @author bl
 */
public class SphericalHarmonicTransform {
    private final Pointer<ComplexDouble> coefficients;
    private final int bandlimit;
    private final boolean reality;
    public Pointer<ComplexDouble> getCoefficients() {
        return coefficients;
    }
    public int getBandlimit() {
        return this.bandlimit;
    }
    public boolean getReality() {
        return this.reality;
    }
    private boolean detectReality() {
        Pointer<Double> pd = coefficients.as(Double.class);
    	for (int el = 0; el < bandlimit; el++) {
            int em = 0;
            int i = el*el + el + em;
            double val = pd.getDoubleAtOffset(i * 2 * 8 + 8);
            if ( val != 0.0 ){
                return(false);
            }
            for (em = 1; em <= el; em++) {
                i = el*el + el + em;
                double val1 = pd.getDoubleAtOffset(i * 2 * 8);
                double val2 = pd.getDoubleAtOffset(i * 2 * 8 + 8);
                int iop = el*el + el - em;
                double val3 = pd.getDoubleAtOffset(iop * 2 * 8);
                double val4 = pd.getDoubleAtOffset(iop * 2 * 8 + 8);
                if ( val3 != Math.pow(-1.0,em) * val1 || val4 != - Math.pow(-1.0,em+1) * val2 ){
                    return(false);
                }
            }
        }
        return(true);
    }
    static public SphericalHarmonicTransform zeros(int bandlimit) {
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        return new SphericalHarmonicTransform(flm);
    }
    static public SphericalHarmonicTransform random(int bandlimit) {
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
	Random randomGenerator = new Random();
	Pointer<Double> pd = flm.as(Double.class);
    	for (int i = 0; i < bandlimit*bandlimit; i++) {
            pd.setDoubleAtOffset(i * 2 * 8, randomGenerator.nextDouble());
            pd.setDoubleAtOffset(i * 2 * 8 + 8, randomGenerator.nextDouble());
        }
        return new SphericalHarmonicTransform(flm);
    }
    static public SphericalHarmonicTransform randomReal(int bandlimit) {
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        Random randomGenerator = new Random();
	Pointer<Double> pd = flm.as(Double.class);
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
        return new SphericalHarmonicTransform(flm);
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
    private SphericalHarmonicTransform(Pointer<ComplexDouble> flm) {
        this.coefficients = flm;
        this.bandlimit = (int)Math.sqrt(flm.getValidElements());
        this.reality = detectReality();
    }
    static public SphericalHarmonicTransform fromHarmonics(Pointer<ComplexDouble> flm) {
        return new SphericalHarmonicTransform(flm);
    }
    static public SphericalHarmonicTransform fromMap(PixelizedMap f) {
        int bandlimit = f.getResolution();
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        f.computeHarmonics(flm);
        return new SphericalHarmonicTransform(flm);
    }   
    static public SphericalHarmonicTransform fromMap(PixelizedMap f, int bandlimit) {
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        f.computeHarmonics(flm);
        return new SphericalHarmonicTransform(flm);
    }
    static public SphericalHarmonicTransform fromWavelets(AxisymmetricWaveletTransform f_wav) {     
        int bandlimit = f_wav.getBandlimit();
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        f_wav.reconstructHarmonics(flm);
        return new SphericalHarmonicTransform(flm);
    }
    
}
