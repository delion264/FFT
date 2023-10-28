#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <chrono>
#include <algorithm>
#include "fft.h"

using namespace std;

class FFT {
public:
  size_t fft_length;
  complex<double> I;

  FFT(size_t len)
  {
    set_length(len);
    I = -1;
    I = sqrt(I);
  }

  ~FFT()
  {
  }

  size_t get_length(void)
  {
    return fft_length;
  }

  void set_length(size_t l)
  {
    fft_length = l;
  }

  vector<complex<double>> zero_pad_prefix(vector<complex<double>> &v) {
    size_t n = v.size();
    vector<complex<double>>::iterator v_start;
    v_start = v.begin();
    v.insert(v_start, n, 0);
    return v;
  }

  vector<complex<double>> zero_pad_suffix(vector<complex<double>> &v) {
    size_t n = v.size();
    vector<complex<double>>::iterator v_end;
    v_end = v.end();
    v.insert(v_end, n, 0);
    return v;
  }

  vector<complex<double>> fft(vector<complex<double>> &input_signal, FFT_DIRECTION direction)
  {

    vector<complex<double>> even_input_idx;
    vector<complex<double>> odd_input_idx;
    vector<complex<double>> even_fft;
    vector<complex<double>> odd_fft;
    size_t n = input_signal.size();
    complex<double> n_root = (double)direction * exp(2 * M_PI * I / (double)n);
    complex<double> omega = 1;
    complex<double> t_k; /* twiddle factor */

    /* Base case */
    if (n == 1)
    {
      return input_signal;
    }

    for (int i = 0; i < (int)n; ++i)
    {
      if ((i % 2) == 0)
      {
        even_input_idx.push_back(input_signal[i]);
      }
      else
      {
        odd_input_idx.push_back(input_signal[i]);
      }
    }

    vector<complex<double>> output_fft(n, 0);

    /* Recursive call */
    even_fft = fft(even_input_idx, direction);
    odd_fft = fft(odd_input_idx, direction);

    /* Stitch together sub-solutions */
    /* DOES THIS COUNTER BIT REVERSAL?
    for (int k = n-1; k > (int)n / 2; --k)
    {
      t_k = omega * odd_fft[k];
      output_fft[k] = even_fft[k] + t_k;
      output_fft[k - (n / 2)] = even_fft[k] - t_k;
      omega = omega * n_root;
    }
    */

    for (int k = 0; k < (int)n / 2; ++k)
    {
      t_k = omega * odd_fft[k];
      output_fft[k] = even_fft[k] + t_k;
      output_fft[k + (n / 2)] = even_fft[k] - t_k;
      omega = omega * n_root;
    }
    reverse(output_fft.begin(), output_fft.end());
    return output_fft;
  }

  vector<complex<double>> cross_corr(vector<complex<double>> &input_a, vector<complex<double>> &input_b) {
    vector<complex<double>> fft_coeff(2*input_a.size(), 0);
    vector<complex<double>>::iterator it = fft_coeff.begin();
    vector<complex<double>> xcorr;
    vector<complex<double>> pad_input_a = zero_pad_prefix(input_a);
    vector<complex<double>> pad_input_b = zero_pad_prefix(input_b);
    vector<complex<double>> fft_a = fft(pad_input_a, FFT_FORWARD);
    vector<complex<double>> fft_b = fft(pad_input_b, FFT_FORWARD);

    for(int i = 0; i < (int)input_a.size(); ++i) {
      *it = fft_a[i]*conj(fft_b[i]);
      ++it;
      /* TO DO: Implement control measure to ensure iterator does not step out of fft_coeff memory */
    }

    xcorr = fft(fft_coeff, FFT_REVERSE);
    return xcorr;
  }
};

int main()
{
  int f_s = 2000000;
  int f_c = 5000;
  complex<double> I = -1;
  I = sqrt(I);
  size_t fft_len = 16384*2;
  vector<double> t;
  vector<complex<double>> signal;
  vector<complex<double>> fft_out;
  vector<double> fft_mag;
  vector<double>::iterator peak_it;
  vector<double> freq_bins;
  int peak_idx;
  double x_k;
  double bin_width = (double)f_s/fft_len;

  for (auto k = 0; k < (int)fft_len; k++)
  {
    t.push_back(k / (double)(f_s));
  }
  for (auto k = 0; k < (int)fft_len; k++)
  {
    signal.push_back(complex<double>(real(exp(2 * M_PI * (double)f_c * I * t[k])), imag(exp(2 * M_PI * (double)f_c * I * t[k]))));
  }

  // auto start_fft_time = chrono::high_resolution_clock::now();
  FFT transform(fft_len);
  fft_out = transform.fft(signal, FFT_FORWARD);

  for (auto k = 0; k < (int)fft_len; k++)
  {
    x_k = abs(fft_out[k]);
    fft_mag.push_back((x_k*x_k)/(fft_len*fft_len));
  }

  // auto end_fft_time = chrono::high_resolution_clock::now();
  // auto fft_duration = chrono::duration_cast<chrono::microseconds>(end_fft_time - start_fft_time);
  // cout << "FFT execution time (us): " << fft_duration.count() << '\n';

  peak_it = max_element(fft_mag.begin(), fft_mag.end());
  peak_idx = peak_it - fft_mag.begin();

  for (auto k = 0; k < (int)fft_len; k++)
  {
    freq_bins.push_back(k * bin_width);
  }

  cout << "FFT peak frequency: " << freq_bins.at(peak_idx) << " Hz / Bin index: " << peak_idx << '\n';
  cout << "Bin width: " << bin_width << " Hz \n";

  return EXIT_SUCCESS;
}
