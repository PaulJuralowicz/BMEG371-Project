#include <immintrin.h>
#include <iostream>
#include <math.h>
#include <chrono>

#define ROW_NUM 21

void simpleNumericalVersion() {
	std::cout << "simple";
	double L = 20.0 / 1000.0;
	double h = 1 * pow(10.0, -3.0);
	double k = 1 * pow(10.0, -8.0);
	int N = (int)(L / h) + 1; //should be 21
	double D = 0.2;
	double Tmax = 0.1;
	int iterations = (int)(Tmax / k) + 1; //should be 10000001
	double w_prev[ROW_NUM] = { 1 };
	double w_curr[ROW_NUM] = { 0 };
	int x = 0;
	double Dkh = (D * k) / (pow(h, 2.0));
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	for (int j = 1; j < iterations; j++) {
		w_curr[0] = w_prev[0] + Dkh * (2 * w_prev[1] - 2 * w_prev[0]);
		for (int i = 1; i < N - 1; i++) {
			w_curr[i] = w_prev[i] + Dkh * (w_prev[i - 1] - 2 * w_prev[i] + w_prev[i + 1]);
		}
		std::copy(std::begin(w_curr), std::end(w_curr), std::begin(w_prev));
	}
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> time_span = t2 - t1;

	std::cout << "It took me " << time_span.count() << " milliseconds.";
	std::cout << std::endl;
}

void intrinNumericalVersion() {
	std::cout << "Intrin";
	double L = 20.0 / 1000.0;
	double h = 1 * pow(10.0, -3.0);
	double k = 1 * pow(10.0, -8.0);
	int N = (int)(L / h) + 1; //should be 21
	double D = 0.2;
	double Tmax = 0.1;
	int iterations = (int)(Tmax / k) + 1; //should be 10000001
	double w_prev[ROW_NUM] = { 1 };
	double w_curr[ROW_NUM] = { 0 };
	int x = 0;
	double Dkh = (D * k) / (pow(h, 2.0));
	__m512d _w_prev, _Dkh, _two, _w_iL1, _w_iG1, _w_sum;
	__m256d _w_prevs, _Dkhs, _twos, _w_iL1s, _w_iG1s, _w_sums;
	_Dkh = _mm512_set1_pd(Dkh);
	_Dkhs = _mm256_set1_pd(Dkh);
	_two = _mm512_set1_pd(-2.0);
	_twos = _mm256_set1_pd(-2.0);
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	for (int j = 1; j < iterations; j++) {
		w_curr[0] = w_prev[0] + Dkh * (2 * w_prev[1] - 2 * w_prev[0]);
		for (int i = 1; i < N-1; i+= 8) {
			_w_prev = _mm512_set_pd(w_prev[i], w_prev[i + 1], w_prev[i + 2], w_prev[i + 3], w_prev[i + 4], w_prev[i + 5], w_prev[i + 6], w_prev[i + 7]);
			_w_iL1 = _mm512_set_pd(w_prev[i - 1], w_prev[i], w_prev[i + 1], w_prev[i + 2], w_prev[i + 3], w_prev[i + 4], w_prev[i + 5], w_prev[i + 6]);
			_w_iG1 = _mm512_set_pd(w_prev[i + 1], w_prev[i + 2], w_prev[i + 3], w_prev[i + 4], w_prev[i + 5], w_prev[i + 6], w_prev[i + 7], w_prev[i+8]);
			_w_sum = _mm512_add_pd(_w_iL1, _w_iG1);
			_w_sum = _mm512_add_pd(_w_sum, _mm512_mul_pd(_w_prev, _two));
			_w_prev = _mm512_add_pd(_w_prev, _mm512_mul_pd(_w_sum,_Dkh));
			//Look into order here...
			//Maybe switch to the copy thing, it might be faster...
			w_curr[i] = double(_w_prev.m512d_f64[7]);
			w_curr[i+1] = double(_w_prev.m512d_f64[6]);
			w_curr[i+2] = double(_w_prev.m512d_f64[5]);
			w_curr[i+3] = double(_w_prev.m512d_f64[4]);
			w_curr[i+4] = double(_w_prev.m512d_f64[3]);
			w_curr[i+5] = double(_w_prev.m512d_f64[2]);
			w_curr[i+6] = double(_w_prev.m512d_f64[1]);
			w_curr[i+7] = double(_w_prev.m512d_f64[0]);
		}
		//0 for boundary condition...
		_w_prevs = _mm256_set_pd(w_prev[17], w_prev[18], w_prev[19], 0.0);
		_w_iL1s = _mm256_set_pd(w_prev[16], w_prev[17], w_prev[18], 0.0);
		_w_iG1s = _mm256_set_pd(w_prev[18], w_prev[19], w_prev[20], 0.0);
		_w_sums = _mm256_add_pd(_w_iL1s, _w_iG1s);
		_w_sums = _mm256_add_pd(_w_sums, _mm256_mul_pd(_w_prevs, _twos));
		_w_prevs = _mm256_add_pd(_w_prevs, _mm256_mul_pd(_w_sums, _Dkhs));
		//Ordering!
		//Maybe switch to the copy thing, it might be faster...
		w_curr[17] = double(_w_prevs.m256d_f64[3]);
		w_curr[18] = double(_w_prevs.m256d_f64[2]);
		w_curr[19] = double(_w_prevs.m256d_f64[1]);
		w_curr[20] = double(_w_prevs.m256d_f64[0]);
		std::copy(std::begin(w_curr), std::end(w_curr), std::begin(w_prev));
	}
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> time_span = t2 - t1;

	std::cout << "It took me " << time_span.count() << " milliseconds.\n";
	std::cout << std::endl;
}

int main() {
	simpleNumericalVersion();
	intrinNumericalVersion();
	return 0;
}

