struct WaveletTransform {
	int bandlimit;
	int firstscale;
	int lastscale;
	int waveletparameter;
	boolean multiresolution;
	boolean reality;
	complex double* wavelets_lm;
};

double* getScale_lm(struct WaveletTransform wav, int j);

int getScaleSize_lm(struct WaveletTransform wav, int j);
int getScaleBandlimit(struct WaveletTransform wav, int j);

void getScale_MW(complex double* wav, struct WaveletTransform wav, int j); // following multires field
void getScale_MW(complex double* wav, struct WaveletTransform wav, int j, int L); // specified explicitly
// what about real signals ?

double* getScale_HPX(struct WaveletTransform wav, int j, int nside); // healpix nsisde
