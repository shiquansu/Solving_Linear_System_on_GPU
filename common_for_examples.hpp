#ifndef __COMMON_FOR_EXAMPLES_HPP__
#define __COMMON_FOR_EXAMPLES_HPP__

 


#include <CL/sycl.hpp>

 

#include <complex>
#include <cstddef>
#include <limits>
#include <list>
#include <map>
#include <type_traits>

 

// for mkl_allocator()
#include "mkl.h"

 


// Half precision support
#if defined(__SYCL_COMPILER_VERSION) && (__SYCL_COMPILER_VERSION < 20190912)
std::ostream &operator<<(std::ostream &out, const sycl::half &v)
{
    out << (double) v;
    return out;
}
#endif

 


//
// custom helpers for referring to different sycl devices
//
enum my_sycl_device_types {cpu_device, gpu_device};

 


std::map<my_sycl_device_types, std::string>
sycl_device_names = { {cpu_device, "CPU"},
                      {gpu_device, "GPU"}  };

 

//
// Users can set environment flags like SYCL_DEVICES_{all,cpu,gpu} to specify
// which devices to run on
//
void set_list_of_devices(std::list<my_sycl_device_types> & list_of_devices)
{
#if defined(SYCL_DEVICES_all)
    list_of_devices.push_back(cpu_device);
    list_of_devices.push_back(gpu_device);
#else
#if defined(SYCL_DEVICES_cpu)
    list_of_devices.push_back(cpu_device);
#endif
#if defined(SYCL_DEVICES_gpu)
    list_of_devices.push_back(gpu_device);
#endif
#endif
}

 

//
// sets device using *_selector_v functionality and returns whether it was successful at finding it
//
void get_sycl_device(sycl::device &my_dev, bool & my_dev_is_found, my_sycl_device_types & desired_sycl_device)
{

 

    my_dev_is_found = true;

 

    try {

 

        switch (desired_sycl_device) {
            case cpu_device:
                my_dev = sycl::device(sycl::cpu_selector_v);
                break;
            case gpu_device:
                my_dev = sycl::device(sycl::gpu_selector_v);
                break;
        }

 

    } catch (...) {
        my_dev_is_found = false;
    }

 

}

 

//
// Sets device using *_selector_v functionality
//
void get_sycl_device(sycl::device &my_dev, my_sycl_device_types & desired_sycl_device)
{
    switch (desired_sycl_device) {
        case cpu_device:
            my_dev = sycl::device(sycl::cpu_selector_v);
            break;
        case gpu_device:
            my_dev = sycl::device(sycl::gpu_selector_v);
            break;
    }
}

 

//
// Get list of available devices to run on
//
std::list<my_sycl_device_types> get_list_of_found_devices(bool failDueToMissingDevice)
{
    std::list<my_sycl_device_types> listOfDevices;
    std::list<my_sycl_device_types> listOfFoundDevices;
    set_list_of_devices(listOfDevices);

 

    failDueToMissingDevice = false;
    for (auto &deviceType : listOfDevices) {
        sycl::device myDev;
        bool myDevIsFound = false;
        get_sycl_device(myDev, myDevIsFound, deviceType);
        if (myDevIsFound) {
            listOfFoundDevices.push_back(deviceType);
        } else {
#ifdef FAIL_ON_MISSING_DEVICES
        std::cout << "No " << sycl_device_names[deviceType] << " devices found; Fail on missing devices is enabled.\n";
        failDueToMissingDevice = true;
#else
        std::cout << "No " << sycl_device_names[deviceType] << " devices found; skipping " << sycl_device_names[deviceType] << " tests.\n";
#endif
        }
    }
    return listOfFoundDevices;
}

 


//
// oneMKL Memory allocator class using mkl_malloc/mkl_free for aligned
// memory allocation/deallocation in STL types like std::vector
//
template <typename T, int align>
struct mkl_allocator
{
    typedef T*          pointer;
    typedef const T*    const_pointer;
    typedef void*       void_pointer;
    typedef const void* const_void_pointer;
    typedef T           value_type;
    typedef size_t      size_type;
    typedef ptrdiff_t   difference_type;

 

    template <typename U> struct rebind { typedef mkl_allocator<U,align> other; };

 

    mkl_allocator() noexcept {}
    template <typename U, int align2> mkl_allocator(mkl_allocator<U,align2> &other)  noexcept {}
    template <typename U, int align2> mkl_allocator(mkl_allocator<U,align2> &&other) noexcept {}

 

    T* allocate(size_t n) {
        void *mem = mkl_malloc(n * sizeof(T), align);
        if (!mem) throw std::bad_alloc();

 

        return static_cast<T*>(mem);
    }

 

    void deallocate(T *p, size_t n) noexcept {
        mkl_free(p);
    }

 

    constexpr size_t max_size() const noexcept {
        return std::numeric_limits<size_t>::max() / sizeof(T);
    }

 

    template <typename U, int align2> constexpr bool operator==(const mkl_allocator<U,align2>) const noexcept { return true;  }
    template <typename U, int align2> constexpr bool operator!=(const mkl_allocator<U,align2>) const noexcept { return false; }

 

    typedef std::true_type is_always_equal;
};

 


//
// print a 2x2 block of data from matrix M using the sycl accessor
//
// M = [ M_00, M_01 ...
//     [ M_10, M_11 ...
//     [ ...
//
template <typename T>
void print_2x2_matrix_values(T M, int ldM, std::string M_name)
{

 

    std::cout << std::endl;
    std::cout << "\t\t\t" << M_name << " = [ " << M[0*ldM + 0] << ", " << M[1*ldM + 0]         << ", ...\n";
    std::cout << "\t\t\t    [ "                << M[0*ldM + 1] << ", " << M[1*ldM + 1] << ", ...\n";
    std::cout << "\t\t\t    [ "                << "...\n";
    std::cout << std::endl;

 

}

 


//
// print a 2x1 block of data from vector x using the sycl accessor
//
// x = [ x_00 ]
//     [ x_10 ]
//     [ ...  ]
//
template <typename fp, sycl::access::target target>
void print_2x1_vector_values(sycl::accessor<fp,1,sycl::access::mode::read,target> &x_accessor, int incx, std::string x_name)
{

 

    std::cout << std::endl;
    std::cout << "\t\t\t" << x_name << " = [ " << x_accessor[0] << " ]\n";
    std::cout << "\t\t\t    [ "                << x_accessor[1*incx] << " ]\n";
    std::cout << "\t\t\t    [ "                << "... ]\n";
    std::cout << std::endl;

 

}

 


//
// print a 1x1 block of data from scalar x using the sycl accessor
//
// x = [ x_00 ]
//

 

template <typename fp, sycl::access::target target>
void print_1x1_scalar_values(sycl::accessor<fp,1,sycl::access::mode::read,target> &x_accessor, std::string x_name)
{

 

    std::cout << std::endl;
    std::cout << "\t\t\t" << x_name << " = [ " << x_accessor[0] << " ]\n";
    std::cout << std::endl;

 

}

 


//
// helpers for initializing templated scalar data type values.
//
template <typename fp>
fp set_fp_value(fp arg1, fp arg2 = 0.0)
{
    return arg1;
}

 

template <typename fp>
std::complex<fp> set_fp_value(std::complex<fp> arg1, std::complex<fp> arg2 = 0.0){
    return std::complex<fp>(arg1.real(), arg2.real());
}

 

bool isDoubleSupported(sycl::device my_dev) {
    return my_dev.get_info<sycl::info::device::double_fp_config>().size() != 0;
}

 

#ifndef NO_MATRIX_HELPERS
//
// Matrix helpers.
//
template <typename T> constexpr T inner_dimension(oneapi::mkl::transpose trans, T m, T n)     { return (trans == oneapi::mkl::transpose::nontrans) ? m : n; }
template <typename T> constexpr T outer_dimension(oneapi::mkl::transpose trans, T m, T n)     { return (trans == oneapi::mkl::transpose::nontrans) ? n : m; }
template <typename T> constexpr T matrix_size(oneapi::mkl::transpose trans, T m, T n, T ldm)  { return outer_dimension(trans, m, n) * ldm; }
template <typename T> constexpr T matrix_size(oneapi::mkl::layout layout_val, T nrows, T ncols, T ld) { return (layout_val == oneapi::mkl::layout::R ? nrows : ncols) * ld; }

 

//
// Random initialization of scalar, vector, general matrix and triangular matrix
//
template <typename fp> fp rand_scalar() { return fp(std::rand()) / fp(RAND_MAX) - fp(0.5); }
template <typename fp> std::complex<fp> rand_complex_scalar() { return std::complex<fp>(rand_scalar<fp>(), rand_scalar<fp>()); }
template <> std::complex<float>  rand_scalar() { return rand_complex_scalar<float>(); }
template <> std::complex<double> rand_scalar() { return rand_complex_scalar<double>(); }
template <> sycl::half rand_scalar()           { return sycl::half(std::rand() % 32000) / sycl::half(32000) - sycl::half(0.5); }

 

template <typename fp> fp rand_scalar(int mag) { fp tmp = fp(mag) + fp(std::rand()) / fp(RAND_MAX) - fp(0.5); if (std::rand() % 2) return tmp; else return -tmp; }
template <typename fp> std::complex<fp> rand_complex_scalar(int mag) { return std::complex<fp>(rand_scalar<fp>(mag), rand_scalar<fp>(mag)); }
template <> std::complex<float>  rand_scalar(int mag) { return rand_complex_scalar<float>(mag); }
template <> std::complex<double> rand_scalar(int mag) { return rand_complex_scalar<double>(mag); }

 

template <typename vec> void rand_vector(vec &v, int n, int inc)
{
    using fp = typename vec::value_type;
    int abs_inc = std::abs(inc);

 

    v.resize(n * abs_inc);

 

    for (int i = 0; i < n; i++)
        v[i*abs_inc] = rand_scalar<fp>();
}

 

template <typename fp> void rand_vector(fp *v, int n, int inc)
{
    int abs_inc = std::abs(inc);
    for (int i = 0; i < n; i++)
        v[i*abs_inc] = rand_scalar<fp>();
}

 


template <typename vec> void rand_matrix(vec &M, oneapi::mkl::layout layout_val, int nrows, int ncols, int ld, int offset = 0, bool resize = true)
{
    using fp = typename vec::value_type;

 

    if (resize)
        M.resize(matrix_size(layout_val, nrows, ncols, ld));

 

    if (layout_val == oneapi::mkl::layout::col_major) {
        for (int j = 0; j < ncols; j++)
            for (int i = 0; i < nrows; i++)
                M[offset + i + j * ld] = rand_scalar<fp>();
    } else {
        for (int i = 0; i < nrows; i++)
            for (int j = 0; j < ncols; j++)
                M[offset + j + i * ld] = rand_scalar<fp>();
    }
}

 

template <typename fp> void rand_matrix(fp *M, oneapi::mkl::layout layout_val, int nrows, int ncols, int ld)
{

 

    if (layout_val == oneapi::mkl::layout::col_major) {
        for (int j = 0; j < ncols; j++)
            for (int i = 0; i < nrows; i++)
                M[i + j * ld] = rand_scalar<fp>();
    } else {
        for (int i = 0; i < nrows; i++)
            for (int j = 0; j < ncols; j++)
                M[j + i * ld] = rand_scalar<fp>();
    }
}

 

template <typename vec> void rand_matrix(vec &M, oneapi::mkl::transpose trans, int m, int n, int ld, int offset = 0, bool resize = true)
{
    using fp = typename vec::value_type;

 

    if (resize)
        M.resize(matrix_size(trans, m, n, ld));

 

    if (trans == oneapi::mkl::transpose::nontrans) {
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                M[offset + i + j * ld] = rand_scalar<fp>();
    } else {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M[offset + j + i * ld] = rand_scalar<fp>();
    }
}

 

template <typename fp> void rand_matrix(fp *M, oneapi::mkl::transpose trans, int m, int n, int ld)
{

 

    if (trans == oneapi::mkl::transpose::nontrans) {
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                M[i + j * ld] = rand_scalar<fp>();
    } else {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M[j + i * ld] = rand_scalar<fp>();
    }
}

 

template <typename vec> void rand_pos_def_matrix(vec &M, int n, int ld)
{
    rand_matrix(M, oneapi::mkl::transpose::nontrans, n, n, ld);
    for (int i = 0; i < n; i++)
        M[i + i * ld] += static_cast<typename vec::value_type>(n);
    return;
}

 

template <typename vec> void rand_trsm_matrix(vec &M, oneapi::mkl::transpose trans, int m, int n, int ld)
{
    using fp = typename vec::value_type;

 

    M.resize(matrix_size(trans, m, n, ld));

 

    if (trans == oneapi::mkl::transpose::nontrans) {
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++) {
                if (i == j)
                    M[i + j * ld] = rand_scalar<fp>(10);
                else
                    M[i + j * ld] = rand_scalar<fp>();
            }
    } else {
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                if (i == j)
                    M[j + i * ld] = rand_scalar<fp>(10);
                else
                    M[j + i * ld] = rand_scalar<fp>();
            }
    }
}

 

template <typename fp> void rand_getrfnp_matrix(fp *M, int m, int n, int ld)
{
    rand_matrix(M, oneapi::mkl::transpose::nontrans, m, n, ld);
    for (int i = 0; i < std::min(m,n); i++)
        M[i + i * ld] += i+1;
    return;
}
#endif  // !NO_MATRIX_HELPERS

 

// TODO: Remove it when compiler stops support for get_cl_code. Keep only .code() usage
template <typename T, typename = void>
struct has_member_code_meta : std::false_type {};

 

template <typename T>
struct has_member_code_meta<T, std::void_t<decltype( std::declval<T>().code() )> > : std::true_type {};

 

template <typename T, typename std::enable_if<has_member_code_meta<T>::value>::type* = nullptr >
auto get_error_code (T x) {
    return x.code().value();
};

 

template <typename T, typename std::enable_if<!has_member_code_meta<T>::value>::type* = nullptr >
auto get_error_code (T x) {
    return x.get_cl_code();
};

 

#endif // __COMMON_FOR_EXAMPLES_HPP__
