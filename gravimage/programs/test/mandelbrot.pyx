# mandelcy1.pyx
# cython: profile=True

import cython

@cython.profile(False)
cdef inline int mandel(double real, double imag, int max_iterations=20):
    '''determines if a point is in the Mandelbrot set based on deciding if,
       after a maximum allowed number of iterations, the absolute value of
       the resulting number is greater or equal to 2.'''
    cdef double z_real = 0., z_imag = 0.
    cdef int i

    for i in range(0, max_iterations):
        z_real, z_imag = ( z_real*z_real - z_imag*z_imag + real,
                           2*z_real*z_imag + imag )
        if (z_real*z_real + z_imag*z_imag) >= 4:
            return i
    return -1

def create_fractal( double min_x,
                    double min_y,
                    double pixel_size,
                    int nb_iterations,
                    colours,
                    image):

    cdef int width, height
    cdef int x, y, start_y, end_y
    cdef int nb_colours, current_colour, new_colour
    cdef double real, imag

    nb_colours = len(colours)
    # image is an ndarray of size: w,h,3
    width = image.shape[0]
    height = image.shape[1]

    for x in range(width):
        real = min_x + x*pixel_size
        for y in range(height):
            imag = min_y + y*pixel_size
            colour = mandel(real, imag, nb_iterations)
            image[x, y, :] = colours[ colour%nb_colours ]
