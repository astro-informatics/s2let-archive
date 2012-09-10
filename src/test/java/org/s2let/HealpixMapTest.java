
package org.s2let;

import org.bridj.Pointer;
import org.junit.*;
import static org.junit.Assert.*;

/**
 *
 * @author bl
 */
public class HealpixMapTest {
    
    @Test
    public void test_Healpix_SphericalHarmonicTransform() {
        boolean verbosity = true;
        test_Healpix_SphericalHarmonicTransform(64, 16, verbosity);
        test_Healpix_SphericalHarmonicTransform(256, 128, verbosity);
    }
    public void test_Healpix_SphericalHarmonicTransform(int nside, int bandlimit, boolean verbosity) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.randomReal(bandlimit);
        PixelizedMap map = PixelizedMap.fromHarmonics(flm, SamplingScheme.HEALPIX, nside);
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromMap(map, bandlimit);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if(verbosity){
            System.out.println("MaxAbsErr of HPX-SHA with L=" + bandlimit + " and nside=" +nside + " : " + error);
        }
        assertEquals(error, 0.0, 1e-1);
    }
    
}
