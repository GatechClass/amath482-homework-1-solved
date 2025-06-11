# amath482-homework-1-solved
**TO GET THIS SOLUTION VISIT:** [AMATH482-Homework 1 Solved](https://mantutor.com/product/amath482-homework-1-solved/)


---

**For Custom/Order Solutions:** **Email:** mantutorcodes@gmail.com  

*We deliver quick, professional, and affordable assignment help.*

---

<h2>Description</h2>



<div class="kk-star-ratings kksr-auto kksr-align-center kksr-valign-top" data-payload="{&quot;align&quot;:&quot;center&quot;,&quot;id&quot;:&quot;81893&quot;,&quot;slug&quot;:&quot;default&quot;,&quot;valign&quot;:&quot;top&quot;,&quot;ignore&quot;:&quot;&quot;,&quot;reference&quot;:&quot;auto&quot;,&quot;class&quot;:&quot;&quot;,&quot;count&quot;:&quot;2&quot;,&quot;legendonly&quot;:&quot;&quot;,&quot;readonly&quot;:&quot;&quot;,&quot;score&quot;:&quot;5&quot;,&quot;starsonly&quot;:&quot;&quot;,&quot;best&quot;:&quot;5&quot;,&quot;gap&quot;:&quot;4&quot;,&quot;greet&quot;:&quot;Rate this product&quot;,&quot;legend&quot;:&quot;5\/5 - (2 votes)&quot;,&quot;size&quot;:&quot;24&quot;,&quot;title&quot;:&quot;AMATH482-Homework 1 Solved&quot;,&quot;width&quot;:&quot;138&quot;,&quot;_legend&quot;:&quot;{score}\/{best} - ({count} {votes})&quot;,&quot;font_factor&quot;:&quot;1.25&quot;}">

<div class="kksr-stars">

<div class="kksr-stars-inactive">
            <div class="kksr-star" data-star="1" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="2" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="3" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="4" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="5" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
    </div>

<div class="kksr-stars-active" style="width: 138px;">
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
    </div>
</div>


<div class="kksr-legend" style="font-size: 19.2px;">
            5/5 - (2 votes)    </div>
    </div>
Abstract

In this paper, Fast Fourier Transform(FFT) is used to find the center frequency and location of a moving submarine by analyzing noisy acoustic data. Each recorded time’s frequency spectrum is averaged to reduce the white noise and obtain the center frequency by applying FFT. I obtained the path of the submarine in time dimension after denoise the data with a Gaussian filter.

1 &nbsp; Introduction

Acoustic data detection has been widely used in security and monitoring. A new submarine technology also utilizes acoustic data to emit an unknown acoustic frequency. By applying FFT, the acoustic data is analyzed to determine a submarine’s center frequency through averaging the spectrum. Besides, the acoustic data is denoised by a Gaussian filter to locate its path. Once figuring out the trajectory of the submarine by receiving the acoustic data, P-8 Poseidon subtracking aircrafts would be sent to track submarines.

2 &nbsp; Theoretical Background

Fourier Series:

Any function in the interval 𝑥 ∈ [−𝐿, 𝐿] can be expressed in terms of linear combination of harmonically sinusoidal functions sin(2𝜋𝑓𝑥 + 𝜙), 𝑐𝑜𝑠(2𝜋𝑓𝑥 + 𝜙) of different amplitudes(A), frequencies(f) and phases(𝜙), known as Fourier series. We write the periodic function as the infinite sum of these sines and cosines:

𝑓 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;(1)

”

Where the coefficient 𝑎’, 𝑎#, and 𝑏# are given by the formulae

𝑎’ &nbsp;𝑑𝑥 &nbsp; &nbsp;𝑎# &nbsp;𝑐𝑜𝑠 #* ( 𝑥𝑑𝑥 &nbsp; &nbsp; 𝑏# &nbsp;𝑠𝑖𝑛 &nbsp;#*( 𝑥𝑑𝑥

Using exponential form of the cos(𝑥) = &nbsp;+”#,+$”# and sin(𝑥) = + “#)+$”# in equation

” &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;”

(1) and coefficient formulae, we can express 𝑓(𝑥) in complex exponential form.

1( &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;)(.$*( )0𝑑𝑥, 𝑘 ∈ ℤ

𝑓(𝑥) = A 𝑐$𝑒 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 𝑤ℎ𝑒𝑟𝑒 𝑐$𝑓(𝑥)𝑒

)1(

Fourier Transform:

The Fourier Transform is an integral transform defined over the interval 𝑥 ∈ [−𝐿, 𝐿].

Fourier transform and its inverse are defined as

𝐹G1 𝑒).$0𝑓(𝑥)𝑑𝑥 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; 𝑓(𝑘) = 1 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;G1 𝑒.$0𝐹(𝑥)𝑑𝑥

√2𝜋 )1 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; √2𝜋 )1

From the definition, we set that the transform is over the area 𝑥 ∈ [−∞, ∞] and our special domain is in the definite main 𝑥 ∈ [−𝐿, 𝐿].

Fast Fourier Transform:

The Fast Fourier Transform (FFT) is the algorithm that transform a function of time into a function of frequency. Comparing to Discrete Fourier Transform (DFT) which has a computing complexity 𝑂(𝑁”), FFT has its advantage that it has relatively low operation complexity of 𝑂(𝑁𝑙𝑜𝑔𝑁).

The FFT algorithm performs more efficiently with large data. Also, MATLAB has a built-in function to perform the FFT for data analysis.

Averaging:

The signal data are filled with large amount of white noise which prevent us from figuring out the signals. Since the white noise could be modeled in a normal distribution and added to each Fourier component of the spectrum. We could average over many realizations in frequency space. Therefore, the noise from each realization will cancel out and we could get the central frequency.

Filtering:

There are many ways to filter the signal data in order to denoise and better detect the signal. A spectral filter around the central frequency is used in this project to remove undesired frequencies and most of the white noise. The Gaussian filter is used here: 𝐹(𝑘) = exp &nbsp;(−𝜏(𝑘 − 𝑘’)”)

This is a Gaussian function with 𝜏 determining the width of the filter and the constant 𝑘’ is the target signal.

3 &nbsp;Algorithm Implementation and Development

In order to find the center frequency generated by the submarine and determine the path of the submarine’s movement, four components have been conducted with algorithm implemented.

1. &nbsp; &nbsp; &nbsp;Data Acquisition and Setup

After the acoustic data within 49 time points is loaded, the frequency components have been setup. Since the Fast Fourier Transform assumes 2𝜋 periodic signals and [−𝐿, 𝐿] is the domain of input data, the frequencies are rescaled by ” *. &nbsp; The 3″( D grid [X, Y, Z] in spatial-time discretized domain and the 3-D grid [Kx, Ky, Kz] in frequency domain have been setup.

2. &nbsp; &nbsp; &nbsp;Spectral averaging

Averaging method with 49 realizations has been applied to average the acoustic data. I applied the 3-D Fourier transform (fftn function) to transform the data from spatial domain to frequency domain and rearrange the data by shifting the zero-frequency component to the center of the array(fftshift function). I got the sum of the data from 49 times of the loop and then average the sum to offset noise. I took the absolute value of the averaged data, since the original results in the matrix contain both real and imaginary parts.

3. &nbsp; &nbsp; &nbsp;Filtering and denoising

First, the Gaussian filter is set with replacing the k0 with the center frequency values in the 3-D. A 3*49 matrix is created to record the 3-D locations in the x, y, and zaxis in each time’s measurement. For each time point, I applied the 3-D Fourier transform(fftn function) on the acoustic data and rearrange (fftshift function) it to the center at zero frequency. Then I applied the filter to the signal data in the frequency domain. Inverse 3-D Fourier transform (ifftn function) is used to acquire the signal data in the time domain. Then I found the indexes for peaks and obtained the 3-D locations. Lastly, I used the plot3 function to plot all of the coordinates recorded to figure out the submarine trajectory.

4. &nbsp; &nbsp; &nbsp;Obtaining the Final Locations

By finding the last column of the locations that record all 49 positions of the submarine, I get the final location of the submarine. Since the subtracking aircraft cannot go underwater, I give the final x and y coordinates for the submarine.

4 &nbsp; Computational Results

The center frequency 𝐾2 = [5.3407, −6.9115, 2.1991]

The coordinates of the submarine at each instance of time were obtained and showed in Table 1. The figure 1 shows the trajectory of the submarine.

Figure 1: The trajectory of the submarine over a 24-hour period in half-hour increments.

Time

x

y

z

Time

x

y

z

1

3.125

0

-8.125

26

-2.8125

5.9375

-0.625

2

3.125

0.3125

-7.8125

27

-3.125

5.9375

-0.3125

3

3.125

0.625

-7.5

28

-3.4375

5.9375

0

4

3.125

1.25

-7.1875

29

-4.0625

5.9375

0.3125

5

3.125

1.5625

-6.875

30

-4.375

5.9375

0.625

6

3.125

1.875

-6.5625

31

-4.6875

5.625

0.9375

7

3.125

2.1875

-6.25

32

-5.3125

5.625

1.25

8

3.125

2.5

-5.9375

33

-5.625

5.3125

1.5625

9

3.125

2.8125

-5.625

34

-5.9375

5.3125

1.875

10

2.8125

3.125

-5.3125

35

-5.9375

5

2.1875

11

2.8125

3.4375

-5

36

-6.25

5

2.5

12

2.5

3.75

-4.6875

37

-6.5625

4.6875

2.8125

13

2.1875

4.0625

-4.375

38

-6.5625

4.375

3.125

14

1.875

4.375

-4.0625

39

-6.875

4.0625

3.4375

15

1.875

4.6875

-3.75

40

-6.875

3.75

3.75

16

1.5625

5

-3.4375

41

-6.875

3.4375

4.0625

17

1.25

5

-3.125

42

-6.875

3.4375

4.375

18

0.625

5.3125

-2.8125

43

-6.875

2.8125

4.6875

19

0.3125

5.3125

-2.5

44

-6.5625

2.5

5

20

0

5.625

-2.1875

45

-6.25

2.1875

5

21

-0.625

5.625

-1.875

46

-6.25

1.875

5.625

22

-0.9375

5.9375

-1.875

47

-5.9375

1.5625

5.625

23

-1.25

5.9375

-1.25

48

-5.3125

1.25

5.9375

24

-1.875

5.9375

-1.25

49

– 5

0.9375

6.5625

25

-2.1875

5.9375

-0.9375

Table 1: the x, y and z coordinates at each of 49 time points over a 24-hour span

In order to send my P-8 Poseidon subtracking aircraft, the last location of the submarine traveling over the 24-hour period is needed to obtain. The location at the 49th instance was found at

[x, y] = [-5.0, 0.9375]

5 &nbsp; Summary and Conclusions

To locate the movement of submarine, we use a broad spectrum recording of acoustics to keep track of the location of submarine in 49 time points. Applying the averaging method to average the raw acoustic data which contain white noise, I find the center frequency, which is the max frequency of where the wave hits the submarine. Apart from averaging the data, I used a Gaussian filter to denoise the data, since it is not necessary to assume the mean of noise throughout time is zero within this method. When we know the location of the central frequency, I could locate the trajectory of the moving submarine over a time period.

Depending on different raw data we get, either averaging or various filters could be applied to denoise the noisy data. The denoising methods with the implication of Fourier transform in the data analysis could be used in many natural science research analyses.

Appendix A. MATLAB Functions

1. &nbsp; &nbsp; &nbsp;y = linespace(x1, x2, n) returns a row vector of n evenly spaced points between x1 and x2.

2. &nbsp; &nbsp; &nbsp;[X, Y, Z] = meshgrid(x, y, z) returns 3-D grid coordinates defined by the vectors x, y, and z.

3. &nbsp; &nbsp; &nbsp;B = reshape(A, sz) reshapes A using size vector sz into B

4. &nbsp; &nbsp; &nbsp;Y = fftn(X) returns the multidimensional Fourier transform of an N-D array using a fast Fourier transform algorithm

5. &nbsp; &nbsp; &nbsp;Y = ifftn(X) returns the multidimensional discrete inverse Fourier transform of an N-D array using a fast Fourier transform algorithm

6. &nbsp; &nbsp; &nbsp;Y = fftshift(X) rearranges a Fourier transform X by shifting the zero-frequency component to the center of the array

7. &nbsp; &nbsp; &nbsp;[row, col] = ind2sub(sz, ind) returns the array row and col containing the equivalent row and column subscripts corresponding to the linear indices ind for a matrix of size sz.

Appendix B. Matlab Code

%% Clean workspace clear all; close all; clc

load(‘subdata.mat’) % Imports the data as the 262144×49

(space by time) matrix called subdata

L = 10; % spatial domain n = 64; % Fourier modes x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x; k = (2*pi/(2*L))*[0:(n/2 – 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);

[Kx,Ky,Kz]=meshgrid(ks,ks,ks); % grid in frequency domain

%% problem 1 – center frequency

Utave = zeros(64, 64, 64)

for j = 1:49

Un(:,:,:)=reshape(subdata(:,j),n,n,n);

Utn = fftshift(fftn(Un)); &nbsp; &nbsp;Utave = Utave + Utn; end

Utave = abs(Utave./49);

[maxVals, indices1] = max(Utave(:));

[kx0, ky0, kz0] = ind2sub([n,n,n], indices1);

Kx0 = Kx(kx0, ky0, kz0)

Ky0 = Ky(kx0, ky0, kz0)

Kz0 = Kz(kx0, ky0, kz0)

%% problem 2 – filter locations = zeros(3, 49); tau = 0.2; filter = exp(-tau*((Kx-Kx0).^2+ (Ky-Ky0).^2+ (Kz-

Kz0).^2)); % Define the filter for j = 1:49

Un(:,:,:)=reshape(subdata(:,j),n,n,n);

Utn = fftshift(fftn(Un));

Unft = filter.*Utn; % Apply the filter to the signal in frequency space

Unf = ifftn(Unft); % the signal in time space

M = max(abs(Unf),[],’all’); &nbsp; &nbsp; indices2 = find(abs(Unf)==M); &nbsp; &nbsp; [a, b, c] = ind2sub([n, n, n], indices2); &nbsp; &nbsp; locations(1, j) = X(a,b,c); &nbsp; &nbsp; locations(2, j) = Y(a,b,c); &nbsp; &nbsp; locations(3, j) = Z(a,b,c); end figure(1) plot3(locations(1, :), locations(2, :), locations(3,:),

‘ko-‘, ‘Linewidth’, 1)

title(‘Submarine Movement trajectory’, ‘Fontsize’, 15) xlabel(‘x’) ylabel(‘y’) zlabel(‘z’) axis([-10 10 -10 10 -10 10]), grid on x49 = locations(1, 49); y49 = locations(2, 49); z49 = locations(3, 49);

L = sprintf(‘The 49th location of the submarine is at %s %d %f.’, x49, y49, z49);

%% Problem 3 – get the table of the 49 2-D positions in each time point num = [1:49] l = [locations(1,:);locations(2,:)]

XY_locations = table(num, l) writetable(XY_locations, ‘xy.csv’) type ‘xy.csv’
