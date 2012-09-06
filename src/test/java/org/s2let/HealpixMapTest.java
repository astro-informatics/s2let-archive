
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
        boolean verbosity = false;
        test_Healpix_SphericalHarmonicTransform(64, 16, verbosity);
        test_Healpix_SphericalHarmonicTransform(256, 128, verbosity);
    }
    public void test_Healpix_SphericalHarmonicTransform(int nside, int bandlimit, boolean verbosity) {
        SphericalHarmonicTransform flm = new SphericalHarmonicTransform(bandlimit);
        flm.randomCoefficiensReal();
        
        HealpixRealMap map = new HealpixRealMap(flm, nside);
        SphericalHarmonicTransform flm_rec = new SphericalHarmonicTransform(map, bandlimit);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if(verbosity){
            System.out.println("MaxAbsErr of HPX-SHA with L=" + bandlimit + " and nside=" +nside + " : " + error);
        }
        assertEquals(error, 0.0, 1e-1);
    }
    
}
