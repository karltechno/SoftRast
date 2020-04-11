# SoftRast
Hobby real-time software rasterizer, requires a CPU with AVX2/FMA.

![Crytek Sponza](https://i.imgur.com/CE8z2bV.jpg)

The goal of the project is to learn more about the graphics pipeline, software rendering and SIMD programming. 

## Features:
- SIMD (AVX2) rasterization/shading.
- Perspective correct interpolation of attributes.
- Fixed point rasterization with 8 bits of sub pixel precision.
- Texture sampling with billinear interpolation and tiled/morton order textures.
- Multithreaded geometry processing and rasterization.
- Sort middle architecture.
- Reverse Z depth buffer (compile time toggleable).
- Mip mapping using screen space partial derivatives.
- No runtime memory allocation (all allocations go through thread local linear allocators with a large upfront allocation).
- Simple OBJ model loader. Creates a binary out of the model and textures for instantaneous loading after the first run.

Various improvements are in todo.txt.

# Assorted references
- [Some notes on Graphics Hardware](http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/gh20121127.pdf) 
- [A Parallel Algorithm for Polygon Rasterization](https://www.cs.drexel.edu/~david/Classes/Papers/comp175-06-pineda.pdf) 
- [A Sorting Classification of Parallel Rendering](http://www.cs.cmu.edu/afs/cs/academic/class/15869-f11/www/readings/molnar94_sorting.pdf)
- [A Modern Approach to Software Rasterization](https://www.semanticscholar.org/paper/A-Modern-Approach-to-Software-Rasterization-Taylor/082e452fea9591687e8c5eea51e6c7c5ff1f18aa)
- [Rasterization on Larabee](http://www.cs.cmu.edu/afs/cs/academic/class/15869-f11/www/readings/abrash09_lrbrast.pdf)
- [Advanced Rasterization](https://web.archive.org/web/20120625103536/http://devmaster.net/forums/topic/1145-advanced-rasterization/)
- [Texture tiling and swizzling](https://fgiesen.wordpress.com/2011/01/17/texture-tiling-and-swizzling/)
- [A trip through the Graphics Pipeline](https://fgiesen.wordpress.com/2011/07/09/a-trip-through-the-graphics-pipeline-2011-index/)
- [Optimizing Software Occlusion Culling](https://fgiesen.wordpress.com/2013/02/17/optimizing-sw-occlusion-culling-index/)
- [High-Performance Software Rasterization on GPUs](https://research.nvidia.com/publication/high-performance-software-rasterization-gpus)
- [Rasterization: a Practical Implementation](https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/overview-rasterization-algorithm)
- [MIP-map Level Selection for Texture Mapping](https://pdfs.semanticscholar.org/9eb3/0298a11e40adf3e3e2ca30b67c3cb90565cf.pdf)
