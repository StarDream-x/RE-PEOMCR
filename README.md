# RE: Perceptual error optimization for Monte Carlo rendering, SIGGRAPH2022

The goal of this repo is to reproduce an interesting paper, 
Perceptual error optimization for Monte Carlo rendering, 
whose main contribution is using halftoning literature and high-frequency spectral
range to improve the outcome of MC-rendering.

Here is the paper [page](https://sampling.mpi-inf.mpg.de/2022-chizhov-perception.html).

## Environment

* gcc 7.5.0
* OpenCV

The same as Games101.

## IntelOpenDenoise

1. cd intelOpenDenoise
2. cd bin
3. convert input.ppm -endian LSB input.pfm
4. /oidnDenoise --hdr input.pfm -o denoise.pfm
5. convert denoise.pfm denoise.jpg
6. convert denoise.jpg denoise.ppm

## Getting Start

1. Build the project if necessary. Notice here is a pre-built version with spp=4, RayTracing.
2. Run RayTracing with options:
   * ./RayTracing 0: to get the traditional average MC-rendering result.
   * ./RayTracing 1 input.ppm: to get the spectrum of input image.
   * ./RayTracing 2 input1.ppm input2.ppm: to get the spectrum of different image.
   * ./RayTracing 3 input.ppm: to get the guassianBlur image of input image with kernel size=3*3 and sigma=sqrt(2/pi).
   * ./RayTracing 4 input.ppm: to get the tiled spectrum image with tile size=32
   * ./RayTracing 5 GT.ppm: to get the optimized image with Ground Truth.

## More

Cornell box is the only model in my project, and more models are expected to
 take into consideration in the near future.