
package org.flaglets.s2let;
import static org.flaglets.s2let.bindings.S2letLibrary.*;
import org.bridj.*;


/**
 *
 * @author bl
 */
public class AxisymmetricWaveletTransform {
    
    protected SphericalHarmonicTransform[] wavelets_lm;
    protected PixelizedMap[] wavelets;
    protected final int waveletParameter, firstScale, scaleCount, bandlimit;
    protected final boolean reality, multiresolution; 
    public boolean getMultiresolution() {
        return multiresolution;
    }
    public int getBandlimit() {
        return bandlimit;
    }
    public int getFirstScale() {
        return firstScale;
    }
    public int getWaveletParameter() {
        return waveletParameter;
    }
    public int getScaleCount() {
        return scaleCount;
    }
    public boolean getReality() {
        return this.reality;
    }
    
    private int getMultiresTotalSize_lm() {
        int total = 0;
        for(int j = firstScale; j <= scaleCount; j++){
            int tempBandlimit = Math.min(getScaleBandlimit(j), bandlimit);
            total += tempBandlimit * tempBandlimit;
	}
        return total;
    }
    private int getMultiresTotalSize() {
        int total = 0;
        for(int j = firstScale; j <= scaleCount; j++){
            int tempBandlimit = Math.min(getScaleBandlimit(j), bandlimit);
            total += tempBandlimit * (2 * tempBandlimit - 1);
	}
        return total;
    }
    public int getScaleBandlimit(int j) {
        return s2let_bandlimit(waveletParameter, j);
    }
    
    public AxisymmetricWaveletTransform(SphericalHarmonicTransform flm, int waveletParameter, int firstScale) {
        this(flm, waveletParameter, firstScale, false);
    }
    public AxisymmetricWaveletTransform(SphericalHarmonicTransform flm, int waveletParameter, int firstScale, boolean multires) {
        this.reality = flm.getReality();
        this.bandlimit = flm.getBandlimit();
        this.scaleCount = s2let_j_max(bandlimit, waveletParameter);
        this.firstScale = firstScale;
        this.waveletParameter = waveletParameter;
        this.multiresolution = multires;
    	Pointer<Double> wav_lm = Pointer.allocateArray(Double.class, (scaleCount + 1) * bandlimit);
    	Pointer<Double> scal_lm = Pointer.allocateArray(Double.class, bandlimit);
    	s2let_axisym_wav_lm(wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        if(!multires){
            Pointer<ComplexDouble> f_wav_lm = Pointer.allocateArray(ComplexDouble.class, (scaleCount + 1 - firstScale) * bandlimit * bandlimit);
            Pointer<ComplexDouble> f_scal_lm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);   
            s2let_axisym_wav_analysis_multires_lm(f_wav_lm, f_scal_lm, flm.getCoefficients(), wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }else{
            Pointer<ComplexDouble> f_wav_lm = Pointer.allocateArray(ComplexDouble.class, getMultiresTotalSize_lm());
            Pointer<ComplexDouble> f_scal_lm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
            s2let_axisym_wav_analysis_lm(f_wav_lm, f_scal_lm, flm.getCoefficients(), wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }
        // EXTRACT INDIVIDUAL WAVELETS IN HARMONIC SPACE
    }
    
    public AxisymmetricWaveletTransform(MWCmplxMap f, int waveletParameter, int firstScale) {
        this(f, waveletParameter, firstScale, false);
    }
    public AxisymmetricWaveletTransform(MWCmplxMap f, int waveletParameter, int firstScale, boolean multires) {
        this.reality = false;
        this.bandlimit = f.getResolution();
        this.scaleCount = s2let_j_max(bandlimit, waveletParameter);
        this.firstScale = firstScale;
        this.waveletParameter = waveletParameter;
        this.multiresolution = multires;      
        if(!multires){
            Pointer<ComplexDouble> f_wav = Pointer.allocateArray(ComplexDouble.class, (scaleCount + 1 - firstScale) * bandlimit * (2 * bandlimit - 1));
            Pointer<ComplexDouble> f_scal = Pointer.allocateArray(ComplexDouble.class, bandlimit * (2 * bandlimit - 1));   
            s2let_axisym_wav_analysis_multires(f_wav, f_scal, f.getMap(), waveletParameter, bandlimit, firstScale);
        }else{
            Pointer<ComplexDouble> f_wav = Pointer.allocateArray(ComplexDouble.class, getMultiresTotalSize());
            Pointer<ComplexDouble> f_scal = Pointer.allocateArray(ComplexDouble.class, bandlimit * (2 * bandlimit - 1));
            s2let_axisym_wav_analysis(f_wav, f_scal, f.getMap(), waveletParameter, bandlimit, firstScale);
        }
        // EXTRACT INDIVIDUAL WAVELETS IN REAL SPACE
    }
    public AxisymmetricWaveletTransform(MWRealMap f, int waveletParameter, int firstScale) {
        this(f, waveletParameter, firstScale, false);
    }  
    public AxisymmetricWaveletTransform(MWRealMap f, int waveletParameter, int firstScale, boolean multires) {
        this.reality = true;
        this.bandlimit = f.getResolution();
        this.scaleCount = s2let_j_max(bandlimit, waveletParameter);
        this.firstScale = firstScale;
        this.waveletParameter = waveletParameter;
        this.multiresolution = multires;        
        if(!multires){
            Pointer<Double> f_wav = Pointer.allocateArray(Double.class, (scaleCount + 1 - firstScale) * bandlimit * (2 * bandlimit - 1));
            Pointer<Double> f_scal = Pointer.allocateArray(Double.class, bandlimit * (2 * bandlimit - 1));   
            s2let_axisym_wav_analysis_multires_real(f_wav, f_scal, f.getMap(), waveletParameter, bandlimit, firstScale);
        }else{
            Pointer<Double> f_wav = Pointer.allocateArray(Double.class, getMultiresTotalSize());
            Pointer<Double> f_scal = Pointer.allocateArray(Double.class, bandlimit * (2 * bandlimit - 1));
            s2let_axisym_wav_analysis_real(f_wav, f_scal, f.getMap(), waveletParameter, bandlimit, firstScale);
        }
        // EXTRACT INDIVIDUAL WAVELETS IN REAL SPACE
    }
    public AxisymmetricWaveletTransform(HealpixRealMap f, int bandlimit, int waveletParameter, int firstScale, boolean multires) {
        this.reality = true;
        this.bandlimit = bandlimit;
        this.scaleCount = s2let_j_max(bandlimit, waveletParameter);
        this.firstScale = firstScale;
        this.waveletParameter = waveletParameter;
        this.multiresolution = multires;
        int nside = f.getResolution();
        Pointer<Double> f_wav = Pointer.allocateArray(Double.class, (scaleCount + 1 - firstScale) * 12 * nside * nside);
        Pointer<Double> f_scal = Pointer.allocateArray(Double.class, 12 * nside * nside);   
        s2let_axisym_hpx_wav_analysis_real(f_wav, f_scal, f.getMap(), nside, waveletParameter, bandlimit, firstScale);
        // EXTRACT INDIVIDUAL WAVELETS IN REAL SPACE
    }
    public AxisymmetricWaveletTransform(HealpixRealMap f, int bandlimit, int waveletParameter, int firstScale) {
        this(f, bandlimit, waveletParameter, firstScale, false);
    }
    
    public SphericalHarmonicTransform getScaleHarm(int scale) {
        if (wavelets_lm == null) {
            computeHarms();
        }
        return wavelets_lm[scale]; // TODO
    }
    public PixelizedMap getScaleMap(int scale) {
        if (wavelets == null) {
            computeMaps();
        }
        return wavelets[scale]; // TODO
    }
    protected void computeHarms() {
        // TODO
    }
    protected void computeMaps() {
        // TODO
    }
    
}