set(lib_src_cpu
    characteristics.c    
    presto_funcs.c
    fresnl.c
    ipmpar.c
    median.c
    dcdflib.c
)

set(lib_src_cuda
    device_main.cu
    host_acceleration.cu
    host_allocate_memory.cu
    host_analysis.cu
    host_debug.cu
    host_get_file_data.cu
    host_get_recorded_data.cu
    host_get_user_input.cu
    host_help.cu
    host_periods.cu
    host_main_function.h
    host_rfi.cu
    host_statistics.cu
    host_stratagy.cu
    host_write_file.cu
    fdas_device.cu
    fdas_host.cu
    fdas_util.cu     
)
