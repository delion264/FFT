#ifndef _FTBP_
#define _FTBP_

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <chrono>
#include <algorithm>
#include <cstddef>

using namespace std;

#define _USE_MATH_DEFINES
#define SAMPLES_PER_BUFFER 1024

enum FFT_DIRECTION
{
  FFT_FORWARD = 1,
  FFT_REVERSE = -1
};

class FTBP {
public:
  size_t fft_length;
  complex<double> I;

  FTBP(size_t len);

  ~FTBP();

  size_t get_length(void);

  void set_length(size_t l);

  vector<complex<double>> zero_pad_prefix(vector<complex<double>> &v);

  vector<complex<double>> zero_pad_suffix(vector<complex<double>> &v);

  vector<complex<double>> compute_dbfs(vector<complex<int16_t>> &v);

  vector<complex<double>> fft(vector<complex<double>> &input_signal, FFT_DIRECTION direction);

  vector<complex<double>> cross_corr(vector<complex<double>> &input_a, vector<complex<double>> &input_b);
};
#endif
