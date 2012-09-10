
package org.s2let;

import org.bridj.*;
import static org.s2let.bindings.S2letLibrary.*;

/**
 *
 * @author bl
 */

public abstract class PixelizedMap {
    protected final SamplingScheme scheme;
    protected final int resolution;
    protected final boolean reality;
    public SamplingScheme getScheme() {
        return scheme;
    }
    public boolean getReality() {
        return this.reality;
    };
    public int getResolution() {
        return this.resolution;
    };
    protected PixelizedMap(SamplingScheme sc, int bandlimit, boolean rea) {
        this.scheme = sc;
        this.reality = rea;
        this.resolution = bandlimit;
    }
    static public PixelizedMap fromHarmonics(SphericalHarmonicTransform flm, SamplingScheme sc, int nside) {
        if (sc == SamplingScheme.MW){        
            if (flm.getReality()) {
                return new MWRealMap(flm);
            } else {
                return new MWCmplxMap(flm);
            }
        } else {
            return new HealpixRealMap(flm, nside);
        }
            
    };
    static public PixelizedMap fromHarmonics(SphericalHarmonicTransform flm, SamplingScheme sc) {
        int defaultNside = flm.getBandlimit();
        return fromHarmonics(flm, sc, defaultNside);
    }
    static public PixelizedMap fromWavelets(AxisymmetricWaveletTransform wavelets, SamplingScheme sc, int nside) {
        if (sc == SamplingScheme.MW){
            if (wavelets.getReality()) {
                return new MWRealMap(wavelets);
            } else {
                return new MWCmplxMap(wavelets);
            }
        } else {
            return new HealpixRealMap(wavelets, nside);
        }
    };
     static public PixelizedMap fromWavelets(AxisymmetricWaveletTransform wavelets, SamplingScheme sc) {
        int defaultNside = wavelets.getBandlimit();
        return fromWavelets(wavelets, sc, defaultNside);
    }
    public abstract void computeHarmonics(Pointer<ComplexDouble> flm);
};


class MWCmplxMap extends PixelizedMap {
    private Pointer<ComplexDouble> mapValues; 
    public Pointer<ComplexDouble> getMap() {
        return mapValues;
    }
    protected MWCmplxMap(SphericalHarmonicTransform flm) {
        super(SamplingScheme.MW, flm.getBandlimit(), false);
        this.mapValues = Pointer.allocateArray(ComplexDouble.class, resolution * (2*resolution-1));
        s2let_mw_alm2map(this.mapValues, flm.getCoefficients(), resolution);
    }
    protected MWCmplxMap(AxisymmetricWaveletTransform wavelets) {
        super(SamplingScheme.MW, wavelets.getBandlimit(), false);
        this.mapValues = Pointer.allocateArray(ComplexDouble.class, resolution * (2*resolution-1));
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, resolution * resolution);
        wavelets.reconstructHarmonics(flm);
        s2let_mw_alm2map(this.mapValues, flm, resolution);
    }
    public void computeHarmonics(Pointer<ComplexDouble> flm){
        s2let_mw_map2alm(flm, this.mapValues, resolution);
    };
}

class MWRealMap extends PixelizedMap {
    private Pointer<Double> mapValues; 
    public Pointer<Double> getMap() {
        return mapValues;
    }
    protected MWRealMap(SphericalHarmonicTransform flm) {
        super(SamplingScheme.MW, flm.getBandlimit(), true);
        this.mapValues = Pointer.allocateArray(Double.class, resolution * (2*resolution-1));
        s2let_mw_alm2map_real(this.mapValues, flm.getCoefficients(), resolution);
    }
    protected MWRealMap(AxisymmetricWaveletTransform wavelets) {
        super(SamplingScheme.MW, wavelets.getBandlimit(), true);
        this.mapValues = Pointer.allocateArray(Double.class, resolution * (2*resolution-1));
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, resolution * resolution);
        wavelets.reconstructHarmonics(flm);
        s2let_mw_alm2map_real(this.mapValues, flm, resolution);
        }
    public void computeHarmonics(Pointer<ComplexDouble> flm){
        s2let_mw_map2alm_real(flm, this.mapValues, resolution);
    };

}

class HealpixRealMap extends PixelizedMap {
    private Pointer<Double> mapValues; 
    public Pointer<Double> getMap() {
        return mapValues;
    }
    protected HealpixRealMap(SphericalHarmonicTransform flm, int res) {
        super(SamplingScheme.HEALPIX, res, true);
        this.mapValues = Pointer.allocateArray(Double.class, 12*res*res);
        s2let_hpx_alm2map_real(this.mapValues, flm.getCoefficients(), this.resolution, flm.getBandlimit());
    }
    protected HealpixRealMap(AxisymmetricWaveletTransform wavelets, int res) {
        super(SamplingScheme.HEALPIX, res, true);
        this.mapValues = Pointer.allocateArray(Double.class, 12*res*res);
        int bandlimit = wavelets.getBandlimit();
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        wavelets.reconstructHarmonics(flm);
        s2let_hpx_alm2map_real(this.mapValues, flm, res, bandlimit);
        }
    public void computeHarmonics(Pointer<ComplexDouble> flm){
        int bandlimit = (int)(Math.sqrt(flm.getValidElements()));
        s2let_hpx_map2alm_real(flm, this.mapValues, resolution, bandlimit);
    };
}


