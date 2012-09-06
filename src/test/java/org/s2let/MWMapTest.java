
package org.s2let;

import org.junit.*;
import static org.junit.Assert.*;


/**
 *
 * @author bl
 */
public class MWMapTest {
    
    @Test
    public void test_MW_SphericalHarmonicTransform() {
        boolean verbosity = false;
        test_MW_SphericalHarmonicTransform(16, verbosity);
        test_MW_SphericalHarmonicTransform(64, verbosity);
        test_MW_SphericalHarmonicTransform(128, verbosity);
    }
    public void test_MW_SphericalHarmonicTransform(int bandlimit, boolean verbosity) {
        SphericalHarmonicTransform flm = new SphericalHarmonicTransform(bandlimit);
        flm.randomCoefficiens();
        MWCmplxMap map = new MWCmplxMap(flm);
        SphericalHarmonicTransform flm_rec = new SphericalHarmonicTransform(map);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if(verbosity){
            System.out.println("MaxAbsErr of MW-SHA with L=" + bandlimit + " : " + error); 
        }
        assertEquals(error, 0.0, 1e-8);
    }
    
    @Test
    public void test_MWReal_SphericalHarmonicTransform() {
        boolean verbosity = false;
        test_MWReal_SphericalHarmonicTransform(16, verbosity);
        test_MWReal_SphericalHarmonicTransform(64, verbosity);
        test_MWReal_SphericalHarmonicTransform(128, verbosity);
    }
    public void test_MWReal_SphericalHarmonicTransform(int bandlimit, boolean verbosity) {
        SphericalHarmonicTransform flm = new SphericalHarmonicTransform(bandlimit);
        flm.randomCoefficiensReal();
        MWRealMap map = new MWRealMap(flm);
        SphericalHarmonicTransform flm_rec = new SphericalHarmonicTransform(map);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if(verbosity){
            System.out.println("MaxAbsErr of RealMW-SHA with L=" + bandlimit + " : " + error);
        }
        assertEquals(error, 0.0, 1e-8);
    }
}
