# Repetitiveness

A little C++20 tool to compute different measures of repetitiveness given an input file:

* the length $n$ of the input,
* the size $\sigma$​ of the alphabet,
* the zeroth-order entropy $\mathcal{H}_0$,
* the number $r$ of BWT runs,
* the number $z_{78}$ of Lempel-Ziv 78 factors,
* the number $z_{77}$ of Lempel-Ziv 77 factors and
* the substring complexity $\delta$.

```
$ ./repetitiveness /scratch/data/pc/dblp.xml 
loading file ...
computing SA ...
computing CST ...
RESULT file=/scratch/data/pc/dblp.xml n=296135875 sigma=98 h0=5.26206 r=41037553 z78=16205171 z77=9576082 delta=4310457.967742
```

## Usage

### Building

```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Wno-dev ..
make
```

### Requirements

This tool requires the [SDSL ](https://github.com/xxsds/sdsl-lite/)to be installed on your system.

Thanks to SDSL, many computations are supported by semi-external data structures. For performance reasons, however, this tool caches the suffix array in $5n$ bytes of RAM to compute $z_{77}$. Furthermore, $z_{78}$ is computed in RAM and requires $\lceil\!\lceil 17 z_{78} \rceil\!\rceil$ bytes of RAM, the hyperceil operator stemming from capacity doubling in `std::vector`.

### License

```
MIT License

Copyright (c) Patrick Dinklage

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### Credits

The computation of $\delta$ has been borrowed from [regindex/substring-complexity ](https://github.com/regindex/substring-complexity), provided also under the MIT license.