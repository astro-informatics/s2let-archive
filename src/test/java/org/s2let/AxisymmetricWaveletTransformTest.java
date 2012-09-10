/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.s2let;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author bl
 */
public class AxisymmetricWaveletTransformTest {
    
    @Test
    public void test_AxisymmetricWaveletTransformTest_harm() {
        boolean verbosity = true;
        test_AxisymmetricWaveletTransformTest_harm(32, 2, 1, false, false, verbosity);
        test_AxisymmetricWaveletTransformTest_harm(32, 3, 1, true, false, verbosity);
        test_AxisymmetricWaveletTransformTest_harm(64, 2, 2, false, true, verbosity);
        test_AxisymmetricWaveletTransformTest_harm(32, 3, 0, true, true, verbosity);
    }
    public void test_AxisymmetricWaveletTransformTest_harm(int bandlimit, int waveletParameter, int firstScale, boolean multires, boolean reality, boolean verbosity) {
        SphericalHarmonicTransform flm;
        if (!reality) {
            flm = SphericalHarmonicTransform.random(bandlimit);
        } else {
            flm = SphericalHarmonicTransform.randomReal(bandlimit);
        }
        AxisymmetricWaveletTransform wavelets = AxisymmetricWaveletTransform.fromHarmonics(flm, waveletParameter, firstScale, multires); 
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromWavelets(wavelets);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if (verbosity) {
            System.out.println("MaxAbsErr of Harm-AxiWav with L=" + bandlimit + ", B="+waveletParameter+", J_min="+firstScale+", multires="+multires+", reality="+reality+" : " + error);
        }
        assertEquals(error, 0.0, 1e-8);
    }
    
    @Test
    public void test_AxisymmetricWaveletTransformTest_MWMap() {
        boolean verbosity = true;
        test_AxisymmetricWaveletTransformTest_MWMap(32, 2, 1, false, false, verbosity);
        test_AxisymmetricWaveletTransformTest_MWMap(64, 3, 1, true, false, verbosity);
        test_AxisymmetricWaveletTransformTest_MWMap(128, 2, 3, false, true, verbosity);
        test_AxisymmetricWaveletTransformTest_MWMap(32, 2, 0, true, true, verbosity);
    }
    public void test_AxisymmetricWaveletTransformTest_MWMap(int bandlimit, int waveletParameter, int firstScale, boolean multires, boolean reality, boolean verbosity) {
        SphericalHarmonicTransform flm;
        if (!reality) {
            flm = SphericalHarmonicTransform.random(bandlimit);
        } else {
            flm = SphericalHarmonicTransform.randomReal(bandlimit);
        }        
        PixelizedMap map = PixelizedMap.fromHarmonics(flm, SamplingScheme.MW);
        AxisymmetricWaveletTransform wavelets = AxisymmetricWaveletTransform.fromMap(map, waveletParameter, firstScale, multires); 
        PixelizedMap map_rec = PixelizedMap.fromWavelets(wavelets, SamplingScheme.MW);
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromMap(map_rec);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if (verbosity) {
            System.out.println("MaxAbsErr of MW-AxiWav with L=" + bandlimit + ", B="+waveletParameter+", J_min="+firstScale+", multires="+multires+", reality="+reality+" : " + error);
        }
        assertEquals(error, 0.0, 1e-8);
    }
    
    @Test
    public void test_AxisymmetricWaveletTransformTest_HealpixMap() {
        boolean verbosity = true;
        test_AxisymmetricWaveletTransformTest_HealpixMap(32, 2, 1, false, 64, verbosity);
        test_AxisymmetricWaveletTransformTest_HealpixMap(64, 3, 1, true, 128, verbosity);
        test_AxisymmetricWaveletTransformTest_HealpixMap(128, 2, 3, false, 128, verbosity);
        test_AxisymmetricWaveletTransformTest_HealpixMap(32, 2, 0, true, 64, verbosity);
    }
    public void test_AxisymmetricWaveletTransformTest_HealpixMap(int bandlimit, int waveletParameter, int firstScale, boolean multires, int nside, boolean verbosity) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.randomReal(bandlimit);
        PixelizedMap map = PixelizedMap.fromHarmonics(flm, SamplingScheme.HEALPIX, nside);
        AxisymmetricWaveletTransform wavelets = AxisymmetricWaveletTransform.fromMap(map, waveletParameter, firstScale, multires); 
        PixelizedMap map_rec = PixelizedMap.fromWavelets(wavelets, SamplingScheme.HEALPIX, nside);
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromMap(map_rec);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if (verbosity) {
            System.out.println("MaxAbsErr of HPX-AxiWav with L=" + bandlimit + ", B="+waveletParameter+", J_min="+firstScale+", multires="+multires+", nside="+nside+" : " + error);
        }
        assertEquals(error, 0.0, 1e-8);
    }
    
}
