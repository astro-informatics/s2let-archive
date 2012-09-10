
package org.s2let;
import static org.s2let.bindings.S2letLibrary.*;
import org.bridj.*;


/**
 *
 * @author bl
 */
public class AxisymmetricWaveletTransform {
    protected Pointer<ComplexDouble> f_wav_lm, f_scal_lm;
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
    public int getScaleBandlimit(int j) {
        return s2let_bandlimit(waveletParameter, j);
    }
    private AxisymmetricWaveletTransform(SphericalHarmonicTransform flm, int waveletParameter, int firstScale, boolean multires) {
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
            f_wav_lm = Pointer.allocateArray(ComplexDouble.class, (scaleCount + 1 - firstScale) * bandlimit * bandlimit);
            f_scal_lm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);   
            s2let_axisym_wav_analysis_lm(f_wav_lm, f_scal_lm, flm.getCoefficients(), wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }else{
            f_wav_lm = Pointer.allocateArray(ComplexDouble.class, getMultiresTotalSize_lm());
            f_scal_lm = Pointer.allocateArray(ComplexDouble.class, s2let_bandlimit(waveletParameter, firstScale-1) );
            s2let_axisym_wav_analysis_multires_lm(f_wav_lm, f_scal_lm, flm.getCoefficients(), wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }
    }
    static public AxisymmetricWaveletTransform fromMap(PixelizedMap f, int waveletParameter, int firstScale, boolean multires) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.fromMap(f);
        return new AxisymmetricWaveletTransform(flm, waveletParameter, firstScale, multires);
    }
    static public AxisymmetricWaveletTransform fromHarmonics(SphericalHarmonicTransform flm, int waveletParameter, int firstScale, boolean multires) {
        return new AxisymmetricWaveletTransform(flm, waveletParameter, firstScale, multires);
    }
    public SphericalHarmonicTransform getNthScaleHarmonics(int scale) {
        int bl;
        if(multiresolution){
            bl = getScaleBandlimit(scale);
        } else {
            bl = this.bandlimit;
        }
        Pointer<ComplexDouble> flm = f_wav_lm.next(scale-firstScale).validElements(bl * bl);
        return SphericalHarmonicTransform.fromHarmonics(flm);
    }
    public PixelizedMap computeNthScaleMap(int scale, SamplingScheme sc, int nside) {
        return PixelizedMap.fromHarmonics(getNthScaleHarmonics(scale), sc, nside);
    }
    
    public PixelizedMap computeNthScaleMap(int scale, SamplingScheme sc) {
        return PixelizedMap.fromHarmonics(getNthScaleHarmonics(scale), sc);
    }
    
    public void reconstructHarmonics(Pointer<ComplexDouble> flm){
        Pointer<Double> wav_lm = Pointer.allocateArray(Double.class, (scaleCount + 1) * bandlimit);
    	Pointer<Double> scal_lm = Pointer.allocateArray(Double.class, bandlimit);
    	s2let_axisym_wav_lm(wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        if(!multiresolution){
            s2let_axisym_wav_synthesis_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }else{
            s2let_axisym_wav_synthesis_multires_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }
    }
}