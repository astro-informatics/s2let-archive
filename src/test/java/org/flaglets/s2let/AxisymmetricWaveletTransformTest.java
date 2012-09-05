/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.flaglets.s2let;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author bl
 */
public class AxisymmetricWaveletTransformTest {
    
    @Test
    public void test_AxisymmetricWaveletTransformTest_HealpixMap() {
        test_AxisymmetricWaveletTransformTest_HealpixMap(64, 16, 2, 0, false);
        test_AxisymmetricWaveletTransformTest_HealpixMap(256, 128, 3, 2, false);
    }
    public void test_AxisymmetricWaveletTransformTest_HealpixMap(int nside, int bandlimit, int waveletParameter, int firstScale, boolean multires) {
        SphericalHarmonicTransform flm = new SphericalHarmonicTransform(bandlimit);
        flm.randomCoefficiensReal();
        HealpixRealMap map = new HealpixRealMap(flm, nside);
        
        AxisymmetricWaveletTransform wavelets = new AxisymmetricWaveletTransform(map, bandlimit, waveletParameter, firstScale, multires);               
        HealpixRealMap map_rec = new HealpixRealMap(wavelets);
        
        SphericalHarmonicTransform flm_rec = new SphericalHarmonicTransform(map_rec, bandlimit);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        System.out.println("MaxAbsErr of HPX-AxiWav with L=" + bandlimit + " and nside=" +nside + " : " + error);
        assertEquals(error, 0.0, 1e-1);
    }
    
    @Test
    public void test_AxisymmetricWaveletTransformTest_MWMap() {
        System.out.println("- Testing ");
        test_AxisymmetricWaveletTransformTest_MWMap(16, 2, 0, false);
    }
    public void test_AxisymmetricWaveletTransformTest_MWMap(int bandlimit, int waveletParameter, int firstScale, boolean multires) {
        SphericalHarmonicTransform flm = new SphericalHarmonicTransform(bandlimit);
        flm.randomCoefficiens();
        MWCmplxMap map = new MWCmplxMap(flm);
        
        AxisymmetricWaveletTransform wavelets = new AxisymmetricWaveletTransform(map, waveletParameter, firstScale, multires);               
        MWCmplxMap map_rec = new MWCmplxMap(wavelets);
        
        SphericalHarmonicTransform flm_rec = new SphericalHarmonicTransform(map_rec);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        System.out.println("MaxAbsErr of MW-AxiWav with L=" + bandlimit + " : " + error);
        assertEquals(error, 0.0, 1e-8);
    }
    
    // TODO : real MW map
    
    
}
