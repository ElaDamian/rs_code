package gps.acquisition;

import cgra.pe.PETrigonometry;

public class Acquisition {

	// Static input parameters
	private static final int FREQ_STEP_COUNT = 11;
	private static final float FREQ_MOD_MIN = -5000;
	private static final float FREQ_MOD_MAX = 5000;
	private static final float FREQ_MOD_STEP = 1000;
	private static final float FREQ_SAMPLING = 400000;
	private static final float ACQ_THRESHOLD = 0.015f;
	private static final float FPI = (float) Math.PI;

	// Results
	private int dopplerShift = 0;
	private int codeShift = 0;

	// Dynamic input parameters
	private final int sCount;

	// Sample storage for unmodified samples
	private final float[] sample_r;
	private final float[] sample_i;
	private int sampleWriteIndex = 0;

	// Storage for intermediate results
	private final float[] res_r_1;
	private final float[] res_i_1;
	private final float[] res_r_2;
	private final float[] res_i_2;

	// Code storage (not transformed)
	private final float[] code_r;
	private final float[] code_i;
	private int codeWriteIndex = 0;

	// Storage for DFT(code) results
	private final float[] code_r_dft;
	private final float[] code_i_dft;

	public Acquisition(int nrOfSamples) {
		this.sCount = nrOfSamples;

		this.sample_r = new float[this.sCount];
		this.sample_i = new float[this.sCount];
		this.res_r_1 = new float[this.sCount];
		this.res_i_1 = new float[this.sCount];
		this.res_r_2 = new float[this.sCount];
		this.res_i_2 = new float[this.sCount];

		this.code_r = new float[this.sCount];
		this.code_i = new float[this.sCount];
		this.code_r_dft = new float[this.sCount];
		this.code_i_dft = new float[this.sCount];
	}
	
	public void enterSample(float real, float imag) {
		sample_r[sampleWriteIndex] = real;
		sample_i[sampleWriteIndex] = imag;
		++sampleWriteIndex;
	}
	
	public void enterCode(float real, float imag) {
		code_r[codeWriteIndex] = real;
		code_i[codeWriteIndex] = imag;
		++codeWriteIndex;
	}
	
	public boolean startAcquisition() {
		final float lFPI = FPI;

		// --- Prepare the code samples by transforming and complex conjugating them
		// OK
		dft(code_r, code_i, code_r_dft, code_i_dft);

		for (int i = 0; i < sCount; ++i)
			code_i_dft[i] *= -1;

		// --- Calculate Smax on the fly while performing all the other operations
		// OK
		float smax = 0;
		float fdAtMax = 0;
		int tauAtMax = 0;

		// --- Loop to visit all required frequency shifts (fd)
		// OK
		for (int j = 0; j < FREQ_STEP_COUNT; ++j)
		{
			// --- Value of fd for this loop iteration
			// OK
			final float fd = FREQ_MOD_MIN + FREQ_MOD_STEP * j;

			// --- Wipe off the carrier (from the given samples) using the current fd
			// OK
			for (int i = 0; i < sCount; ++i)
			{
				final float angle = 2 * lFPI * fd * i / FREQ_SAMPLING;
				final float sined = PETrigonometry.sin(angle);
				final float cosed = PETrigonometry.cos(angle);

				res_r_1[i] =  sample_r[i] * cosed + sample_i[i] * sined;
				res_i_1[i] = -sample_r[i] * sined + sample_i[i] * cosed;
			}

			// ... Samples in res1

			// --- Apply the DFT to the given samples which had the carrier wiped off
			// OK
			dft(res_r_1, res_i_1, res_r_2, res_i_2);

			// ... Samples in res2

			// --- Multiply both (samples and code) DFT results
			// OK
			for (int i = 0; i < sCount; ++i) {
				res_r_1[i] = res_r_2[i] * code_r_dft[i] - res_i_2[i] * code_i_dft[i];
				res_i_1[i] = res_i_2[i] * code_r_dft[i] + res_r_2[i] * code_i_dft[i];
			}

			// ... Samples in res1

			// --- Apply the IDFT to retrieve the results
			// OK
			idft(res_r_1, res_i_1, res_r_2, res_i_2);

			// ... Samples in res2

			// --- Search the maximum in the resulting vector
			// OK
			for (int i = 0; i < sCount; ++i) {
				final float abs_sqr = res_r_2[i] * res_r_2[i] + res_i_2[i] * res_i_2[i];
				if (smax < abs_sqr) {
					smax = abs_sqr;
					fdAtMax = fd;
					tauAtMax = i;
				}
			}
		}

		// --- Save results
		// OK
		dopplerShift =  (int) fdAtMax;
		codeShift = tauAtMax;

		// --- Calculate the signals' power (do not forget the normalisation)
		// OK
		float pin = 0;
		
		for (int i = 0; i < sCount; ++i)
			pin += sample_r[i] * sample_r[i] + sample_i[i] * sample_i[i];
		
		pin /= sCount;

		// --- Calculate the threshold by using the results
		// OK
		final float resultingThreshold = smax / pin;

		return resultingThreshold > ACQ_THRESHOLD;
	}

	/**
	 * Perform the DFT for the given array of complex float values.
	 * Assumes that all arrays are of equal length N.
	 * Initial version of tho code taken from:
	 * https://www.nayuki.io/page/how-to-implement-the-discrete-fourier-transform
	 *
	 * @param in_r real part of the input values
	 * @param in_i imaginary part of the input values
	 * @param out_r real part of the output values
	 * @param out_i imaginary part of the output values
	 */
	static private void dft(float[] in_r , float[] in_i,
					 		float[] out_r, float[] out_i)
	{
		final float lFPI = FPI;
		final int n = in_r.length;

		for (int k = 0; k < n; k++) { // For each output element
			float sumreal = 0;
			float sumimag = 0;
			for (int t = 0; t < n; t++) { // For each input element
				float angle = 2 * lFPI * t * k / n;
				sumreal +=  in_r[t] * PETrigonometry.cos(angle) + in_i[t] * PETrigonometry.sin(angle);
				sumimag += -in_r[t] * PETrigonometry.sin(angle) + in_i[t] * PETrigonometry.cos(angle);
			}
			out_r[k] = sumreal;
			out_i[k] = sumimag;
		}
	}

	/**
	 * Perform the IDFT for the given array of complex float values.
	 * Assumes that all arrays are of equal length N.
	 *
	 * @param in_r real part of the input values
	 * @param in_i imaginary part of the input values
	 * @param out_r real part of the output values
	 * @param out_i imaginary part of the output values
	 */
	static private void idft(float[] in_r , float[] in_i,
							 float[] out_r, float[] out_i)
	{
		final float lFPI = FPI;
		final int n = in_r.length;

		for (int k = 0; k < n; k++) { // For each output element
			float sumreal = 0;
			float sumimag = 0;
			for (int t = 0; t < n; t++) { // For each input element
				float angle = 2 * lFPI * t * k / n;
				sumreal += in_r[t] * PETrigonometry.cos(angle) - in_i[t] * PETrigonometry.sin(angle);
				sumimag += in_r[t] * PETrigonometry.sin(angle) + in_i[t] * PETrigonometry.cos(angle);
			}
			out_r[k] = sumreal / n;
			out_i[k] = sumimag / n;
		}
	}
	
	public int getDopplerverschiebung(){
		return dopplerShift;
	}
	
	public int getCodeVerschiebung(){
		return codeShift;
	}

}
