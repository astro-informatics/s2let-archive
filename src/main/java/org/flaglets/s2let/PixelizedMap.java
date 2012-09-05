
package org.flaglets.s2let;

import org.bridj.*;
import static org.flaglets.s2let.bindings.S2letLibrary.*;

/**
 *
 * @author bl
 */
public abstract class PixelizedMap {
    protected SamplingScheme scheme;
    protected boolean reality;
    protected int resolution;
    public SamplingScheme getScheme() {
        return scheme;
    }
    public boolean getReality() {
        return this.reality;
    };
    public int getResolution() {
        return this.resolution;
    };
    public PixelizedMap() {
        
    };
}

class MWCmplxMap extends PixelizedMap {
    private Pointer<ComplexDouble> mapValues;
    public Pointer<ComplexDouble> getMap() { 
        return this.mapValues; 
    }
    public long getPixelCount() {
        return this.mapValues.getValidElements();
    };
    public MWCmplxMap(SphericalHarmonicTransform flm) {
        this.resolution = flm.getBandlimit();
        this.scheme = SamplingScheme.MW;
        this.mapValues = Pointer.allocateArray(ComplexDouble.class, resolution * (2*resolution-1)); ;
        s2let_mw_alm2map(this.mapValues, flm.getCoefficients(), resolution);
    }
    public MWCmplxMap(AxisymmetricWaveletTransform wavelets) {
        int bandlimit = wavelets.getBandlimit();
        this.resolution = bandlimit;this.scheme = SamplingScheme.MW;
        this.mapValues = Pointer.allocateArray(ComplexDouble.class, resolution * (2*resolution-1));
        int scaleCount = wavelets.getScaleCount();
        int firstScale = wavelets.getFirstScale();
        boolean multiresolution = wavelets.getMultiresolution();
        int waveletParameter = wavelets.getWaveletParameter();
        if(!multiresolution){
            Pointer<ComplexDouble> f_wav = Pointer.allocateArray(ComplexDouble.class, (scaleCount + 1 - firstScale) * bandlimit * (2 * bandlimit - 1));
            Pointer<ComplexDouble> f_scal = Pointer.allocateArray(ComplexDouble.class, bandlimit * (2 * bandlimit - 1));   
            // Construct right pointer
            s2let_axisym_wav_analysis_multires(mapValues, f_wav, f_scal, waveletParameter, bandlimit, firstScale);
        }else{
            Pointer<ComplexDouble> f_wav = Pointer.allocateArray(ComplexDouble.class, getMultiresTotalSize());
            Pointer<ComplexDouble> f_scal = Pointer.allocateArray(ComplexDouble.class, bandlimit * (2 * bandlimit - 1));
            // Construct right pointer
            s2let_axisym_wav_synthesis(mapValues, f_wav, f_scal, waveletParameter, bandlimit, firstScale);
        }
    }
}

class MWRealMap extends PixelizedMap {
    private Pointer<Double> mapValues; 
    public Pointer<Double> getMap() { 
        return this.mapValues; 
    }
    public long getPixelCount() {
        return this.mapValues.getValidElements();
    };
    public MWRealMap(SphericalHarmonicTransform flm) {
        this.resolution = flm.getBandlimit();
        this.scheme = SamplingScheme.MW;
        this.mapValues = Pointer.allocateArray(Double.class, resolution * (2*resolution-1));
        s2let_mw_alm2map_real(this.mapValues, flm.getCoefficients(), resolution);
    }
    public MWRealMap(AxisymmetricWaveletTransform wavelets) {
        int bandlimit = wavelets.getBandlimit();
        this.resolution = bandlimit;this.scheme = SamplingScheme.MW;
        this.mapValues = Pointer.allocateArray(Double.class, resolution * (2*resolution-1));
        int scaleCount = wavelets.getScaleCount();
        int firstScale = wavelets.getFirstScale();
        boolean multiresolution = wavelets.getMultiresolution();
        int waveletParameter = wavelets.getWaveletParameter();
        if(!multiresolution){
            Pointer<Double> f_wav = Pointer.allocateArray(Double.class, (scaleCount + 1 - firstScale) * bandlimit * (2 * bandlimit - 1));
            Pointer<Double> f_scal = Pointer.allocateArray(Double.class, bandlimit * (2 * bandlimit - 1));   
            // Construct right pointer
            s2let_axisym_wav_analysis_multires_real(mapValues, f_wav, f_scal, waveletParameter, bandlimit, firstScale);
        }else{
            Pointer<Double> f_wav = Pointer.allocateArray(Double.class, getMultiresTotalSize());
            Pointer<Double> f_scal = Pointer.allocateArray(Double.class, bandlimit * (2 * bandlimit - 1));
            // Construct right pointer
            s2let_axisym_wav_synthesis_real(mapValues, f_wav, f_scal, waveletParameter, bandlimit, firstScale);
        }
    }
}

class HealpixRealMap extends PixelizedMap {
    private Pointer<Double> mapValues;
    public Pointer<Double> getMap() { 
        return this.mapValues; 
    }
    public long getPixelCount() {
        return this.mapValues.getValidElements();
    };
    public HealpixRealMap(SphericalHarmonicTransform flm, int res) {
        this.resolution = res;
        this.scheme = SamplingScheme.HEALPIX;
        int bandlimit = flm.getBandlimit();
        this.mapValues = Pointer.allocateArray(Double.class, 12*res*res);
        s2let_hpx_alm2map_real(this.mapValues, flm.getCoefficients(), this.resolution, bandlimit);
    }
    public HealpixRealMap(AxisymmetricWaveletTransform wavelets) {
        // Calculate wavelet decomposition
    }
}

enum SamplingScheme {
    HEALPIX, MW
};

