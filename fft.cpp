#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <complex>
#include "fft.h"

using namespace std;

class FFT {
public:
  size_t fft_length;
  FFT_DIRECTION direction;
  complex<double> I;

  FFT(size_t len, FFT_DIRECTION dir) {
    direction = dir;
    fft_length = len;
    I = -1;
    I = sqrt(I);
  }

  ~FFT() {

  }

  size_t get_length(void) {
    return fft_length;
  }

  FFT_DIRECTION get_direction(void) {
    return direction;
  }

  void set_length(size_t l) {
    fft_length = l;
  }

  void set_direction(FFT_DIRECTION d) {
    direction = d;
  }

  vector<complex<double>> fft(vector<complex<double>>& input_signal) {

    vector<complex<double>> even_input_idx;
    vector<complex<double>> odd_input_idx;
    vector<complex<double>> even_fft;
    vector<complex<double>> odd_fft;

    size_t n = input_signal.size();

    complex<double> n_root = (double)direction*exp(2*M_PI*I/(double)n);
    complex<double> omega = 1;

    /* Base case */
    if (n == 1) {
      return input_signal;
    }

    for(int i = 0; i < n; ++i) {
      if((i % 2) == 0) {
        even_input_idx.push_back(input_signal[i]);
      } else {
        odd_input_idx.push_back(input_signal[i]);
      }
    }

    vector<complex<double>> output_fft (n, 0);

    /* Recursive call */
    even_fft = fft(even_input_idx);
    odd_fft = fft(odd_input_idx);

    /* Stitch together sub-solutions */
    for(int k = 0; k < n/2; ++k) {
      output_fft[k] = even_fft[k] + omega * odd_fft[k];
      output_fft[k + (n/2)] = even_fft[k] - omega * odd_fft[k];
      omega = omega * n_root;
    }

    return output_fft;
  }
};
